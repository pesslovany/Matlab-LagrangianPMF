%%
%   Lagrangian based grid-based filter for real TAN data.
%   authors: Jakub Matousek, pesslovany@gmail.com
%            Jindrich Dunik, dunikj@kky.zcu.cz, University of West Bohemia
% See paper: TBD doi.

%% Parameters and system simulation
clc
clear variables
close all
format shortG

modelChoose = 4; % choose model 3D - 3, 4D with 2D measurement - 4

load('data.mat') % map of terrain
% Map interpolant for measurement
vysky = griddedInterpolant(souradniceX',souradniceY',souradniceZ',"linear","none");
% Time steps
timeSteps = 1:1:length(souradniceGNSS);
% Last time
endTime = length(timeSteps);

% Number of Monte Carlo simulations
MC = 1;

% For cycle over Monte Carlo simulations
for mc = 1:1:MC

    model = initModel(modelChoose, souradniceGNSS, hB, vysky);

    % Unpack model structure variables
    fields = fieldnames(model);
    for i = 1:numel(fields)
        eval([fields{i} ' = model.' fields{i} ';']);
    end

    clear model

    % Example of trajectories where UKF fails
    %load("data_UKFdivergence_3D.mat")
    %load("data_UKFdivergence_4D.mat")
    
    % PF init and param
    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; % initial condition for PF

    % UKF init
    xp_ukf = meanX0;
    Pp_ukf = varX0;
    
    % PMF init and param
    % Initial grid
    [predGrid, predGridDelta, gridDimOld, gridCenter, gridRotation] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    % Initial PMD
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predPdf = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)
    % Auxiliary variables
    predDenDenomW = sqrt((2*pi)^nx*det(Q)); % Denominator for convolution in predictive step
    halfGrid = ceil(N/2); % Middle row of the TPM matrix index
    
    for k = 1:1:endTime-1
        disp(['Step:', num2str(k),'/',num2str(endTime-1)])

        %% Grid-based Filter
        tic
        % Measurement update
        [filtPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predPdf,predGridDelta(:,k),hfunct); % Measurement update           
        % Filtering mean and var
        filtMeanPMF(:,k) = predGrid*filtPdf*prod(predGridDelta(:,k)); % Measurement update mean
        chip_ = (predGrid-filtMeanPMF(:,k));
        chip_w = chip_.*repmat(filtPdf',nx,1);
        filtVarPMF(:,:,k) = chip_w*chip_' * prod(predGridDelta(:,k)); % Measurement update variance
        filtGrid = predGrid;
        % Meas PDF interp
        [filtPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(filtPdf,...
            gridDimOld,gridCenter,gridRotation,Q,sFactor,nx,Npa,filtMeanPMF(:,k),filtVarPMF(:,:,k),F);
        gridCenter = F*filtMeanPMF(:,k) + u(:,k);
        gridRotation = F*(eigVect);
        % Time Update
        [predPdf,predGrid,predGridDelta] = gbfTimeUpdateFFT(F,filtPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));
        tocPMF(k) = toc; % Time evaluation

        %% Partcle Filter
        tic
        %Measurement Update
        predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); % Prediction density grid through measurement EQ
        predThrMeasEq(1,isnan(predThrMeasEq(1,:))) = inf; % To protect from map extrapolations
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        w = pdf(V.pdf,pom); % Weights
        w = w/sum(w); % Normalization of weights
        filtMeanPF(:,k) = ksiPrior*w; % mean filtering estimate
        filtVarPF(:,k) = ksiPrior.^2*w - filtMeanPF(:,k).^2; % diagonal of filtering covariance matrix
        % Resampling
        cumW = cumsum(w);
        randomN = rand(1,noPart);
        % Use this in case the mex file not working
        % for ind = 1:1:noParts
        %     I(ind) = find(cumW >= randomN(ind),1, 'first');
        % end
        I = binarySearch(cumW, randomN, 'first');
        % Time Update
        ksi = ksiPrior(:,I);
        ksiPrior = F*ksi + mvnrnd(zeros(1,nx),Q,noPart)' + u(:,k);
        tocPF(k) = toc;

        % Unscented Kalman Filter (basic/recursive/Joseph version)
        if ukfOn == 1
            tic;
            %% UKF - to compare
            % -- measurement update
            [chip, wm, wc] = msp(xp_ukf(:,k),Pp_ukf(:,:,k),kappa);
            % --
            dzetap = zeros(nz,SPnum);
            for j=1:1:SPnum      
                dzetap(:,j) = hfunct(chip(:,j),zeros(nz,1),k);
            end
            zp = zeros(nz,1);
            for j=1:1:SPnum
              zp = zp + wm(j)*dzetap(:,j);
            end
            zp = zp + meanV;
            % --
            Pzp = zeros(nz);
            for j=1:1:SPnum
              Pzp = Pzp + wc(j)*(dzetap(:,j)-zp)*(dzetap(:,j)-zp)';
            end
            Pep = Pzp + R;    
            % --
            Pxzp = zeros(nx,nz);
            for j=1:1:SPnum
              Pxzp = Pxzp + wc(j)*(chip(:,j)-xp_ukf(:,k))*(dzetap(:,j)-zp)';
            end
            % filtering update
            K = Pxzp*inv(Pep);
            e = z(:,k) - zp;
            xf_ukf(:,k) = xp_ukf(:,k) + K*e;
            Pf_ukf(:,:,k) = Pp_ukf(:,:,k) - Pxzp*inv(Pep)*Pxzp';
            %Happrox = Pxzp'/Pp_ukf(:,:,k);
            %Pf_ukf(:,:,k) = (eye(nx)-K*Happrox)*Pp_ukf(:,:,k)*(eye(nx)-K*Happrox)'+K*R*K';
            % -- time update (Kalman filter)
            xp_ukf(:,k+1) = F*xf_ukf(:,k) + u(:,k);
            Pp_ukf(:,:,k+1) = F*Pf_ukf(:,:,k)*F' + Q;
            tocUKF(k) = toc;
        end

        endtime_act = k;


    end

    % Evaluation data preparation
    rmsePMF(:,mc) = sqrt(mean((x(:,1:k-1)-filtMeanPMF(:,1:k-1)).^2,2)); %#ok<*SAGROW>
    astdPMF11(mc) = sqrt(mean(filtVarPMF(1,1,1:k-1)));
    astdPMF22(mc) = sqrt(mean(filtVarPMF(2,2,1:k-1)));
    astdPMF33(mc) = sqrt(mean(filtVarPMF(3,3,1:k-1)));
    annes_PMF = 0;
    for indAn = 1:1:k-1
        annes_PMF = annes_PMF + ((x(:,indAn)-filtMeanPMF(:,indAn)).*(1./(diag(filtVarPMF(:,:,indAn)))))'*(x(:,indAn)-filtMeanPMF(:,indAn));
    end
    annes_PMFout(mc) = annes_PMF;

    rmsePF(:,mc) = sqrt(mean((x(:,1:k-1)-filtMeanPF(:,1:k-1)).^2,2));
    astdPF11(mc) = sqrt(mean(filtVarPF(1,1:k-1)));
    astdPF22(mc) = sqrt(mean(filtVarPF(2,1:k-1)));
    astdPF33(mc) = sqrt(mean(filtVarPF(3,1:k-1)));
    annes_PF = 0;
    for indAn = 1:1:k-1
        annes_PF = annes_PF + ((x(:,indAn)-filtMeanPF(:,indAn)).*(1./(diag(filtVarPMF(:,:,indAn)))))'*(x(:,indAn)-filtMeanPF(:,indAn));
    end
    annes_PFout(mc) = annes_PF;

    if ukfOn == 1
        rmseUKF(:,mc) = sqrt(mean((x(:,1:k-1)-xf_ukf(:,1:k-1)).^2,2));
        astdUKF11(mc) = sqrt(mean(squeeze(Pf_ukf(1,1,1:k-1))));
        astdUKF22(mc) = sqrt(mean(squeeze(Pf_ukf(2,2,1:k-1))));
        astdUKF33(mc) = sqrt(mean(squeeze(Pf_ukf(3,3,1:k-1))));
        annes_UKF = 0;
        for indAn = 1:1:k-1
            annes_UKF = annes_UKF + ((x(:,indAn)-xf_ukf(:,indAn)).*(1./(diag(Pf_ukf(:,:,indAn)))))'*(x(:,indAn)-xf_ukf(:,indAn));
        end
        annes_UKFout(mc) = annes_UKF;
    end

    tocPMFavg(:,mc) = mean(tocPMF);
    tocPFavg(:,mc) = mean(tocPF);
    tocUKFavg(:,mc) = mean(tocUKF);

end

% Evaluate results, create table and plots
if ukfOn == 0
    evalRes(nx, k, annes_PMFout, annes_PFout, rmsePMF, rmsePF, ...
        astdPMF11, astdPMF22, astdPF11, astdPF22, tocPMFavg, tocPFavg, ...
        x, filtMeanPMF, filtVarPMF, filtMeanPF, filtVarPF)
else
    evalRes_withUKF(nx, k, annes_PMFout, annes_PFout, annes_UKFout, rmsePMF, rmsePF, rmseUKF, ...
    astdPMF11, astdPMF22, astdPF11, astdPF22, astdUKF11, astdUKF22, tocPMFavg, tocPFavg, tocUKFavg, ...
    x, filtMeanPMF, filtVarPMF, filtMeanPF, filtVarPF, xf_ukf, Pf_ukf)
end