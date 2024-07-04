%%
%   Lagrangian based grid-based filter for real TAN data.
%   authors: Jakub Matousek, pesslovany@gmail.com
%            Jindrich Dunik, dunikj@kky.zcu.cz, University of West Bohemia

%% Parameters and system simulation
clc
clear variables
close all
format shortG

modelChoose = 3; % choose model 3, 4D

load('data.mat') % map of terrain
% Map interpolant for measurement
vysky = griddedInterpolant(souradniceX',souradniceY',souradniceZ',"linear","nearest");
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
    
    % PF init and param
    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; % initial condition for PF
    
    % PMF init and param
    % Initial grid
    [predGrid, predGridDelta, gridDimOld, gridCenter, gridRotation] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    % Initial PMD
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)
    % Auxiliary variables
    predDenDenomW = sqrt((2*pi)^nx*det(Q)); % Denominator for convolution in predictive step
    halfGrid = ceil(N/2); % Middle row of the TPM matrix index
    
    for k = 1:1:endTime-1
        disp(['Step:', num2str(k),'/',num2str(endTime-1)])
        %% Grid-based Filter
        tic
        % Measurement update
        [measPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predDensityProb,predGridDelta(:,k),hfunct); % Measurement update
        measPdf(measPdf<1e-11) = 0; %         
        % Filtering mean and var
        measMeanPMF(:,k) = predGrid*measPdf*prod(predGridDelta(:,k)); % Measurement update mean
        chip_ = (predGrid-measMeanPMF(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVarPMF(:,:,k) = chip_w*chip_' * prod(predGridDelta(:,k)); % Measurement update variance
        measGrid = predGrid;
        % Meas PDF interp
        [measPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(measPdf,...
            gridDimOld,gridCenter,gridRotation,Q,sFactor,nx,Npa,measMeanPMF(:,k),measVarPMF(:,:,k),F);
        gridCenter = F*measMeanPMF(:,k) + u(:,k);
        gridRotation = F*(eigVect);
        % Time Update
        [predDensityProb,predGrid,predGridDelta] = gbfTimeUpdateFFT(F,measPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));
        tocPMF(k) = toc; % Time evaluation

        %% Partcle Filter
        tic
        %Measurement Update
        predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); % Prediction density grid through measurement EQ
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        w = pdf(V.pdf,pom); % Weights
        w = w/sum(w); % Normalization of weights
        measMeanPF(:,k) = ksiPrior*w; % mean filtering estimate
        measVarPF(:,k) = ksiPrior.^2*w - measMeanPF(:,k).^2; % diagonal of filtering covariance matrix
        % Resampling
        cumW = cumsum(w);
        randomN = rand(1,noPart);
        I = binarySearch(cumW, randomN, 'first');
        % Time Update
        ksi = ksiPrior(:,I);
        ksiPrior = F*ksi + mvnrnd(zeros(1,nx),Q,noPart)' + u(:,k);
        tocPF(k) = toc;
        kf_act = k;
    end

    % Evaluation data preparation
    rmsePMF(:,mc) = sqrt(mean((x(:,1:k-1)-measMeanPMF(:,1:k-1)).^2,2)); %#ok<*SAGROW>
    astdPMF11(mc) = sqrt(mean(measVarPMF(1,1,1:k-1)));
    astdPMF22(mc) = sqrt(mean(measVarPMF(2,2,1:k-1)));
    astdPMF33(mc) = sqrt(mean(measVarPMF(3,3,1:k-1)));
    annes_PMF = 0;
    for indAn = 1:1:k-1
        annes_PMF = annes_PMF + ((x(:,indAn)-measMeanPMF(:,indAn)).*(1./(diag(measVarPMF(:,:,indAn)))))'*(x(:,indAn)-measMeanPMF(:,indAn));
    end
    annes_PMFout(mc) = annes_PMF;

    rmsePF(:,mc) = sqrt(mean((x(:,1:k-1)-measMeanPF(:,1:k-1)).^2,2));
    astdPF11(mc) = sqrt(mean(measVarPF(1,1:k-1)));
    astdPF22(mc) = sqrt(mean(measVarPF(2,1:k-1)));
    astdPF33(mc) = sqrt(mean(measVarPF(3,1:k-1)));
    annes_PF = 0;
    for indAn = 1:1:k-1
        annes_PF = annes_PF + ((x(:,indAn)-measMeanPF(:,indAn)).*(1./(diag(measVarPMF(:,:,indAn)))))'*(x(:,indAn)-measMeanPF(:,indAn));
    end
    annes_PFout(mc) = annes_PF;

    tocPMFavg(:,mc) = mean(tocPMF);
    tocPFavg(:,mc) = mean(tocPF);

end

% Evaluate results, create table and plots
evalRes(nx, k, annes_PMFout, annes_PFout, rmsePMF, rmsePF, ...
    astdPMF11, astdPMF22, astdPF11, astdPF22, tocPMFavg, tocPFavg, ...
    x, measMeanPMF, measVarPMF, measMeanPF, measVarPF)
