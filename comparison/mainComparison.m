%%
%   Lagrangian based grid-based filter for real TAN data.
%   Comparison with state of the art filters
%   authors: Jakub Matousek, pesslovany@gmail.com
%            Jindrich Dunik, dunikj@kky.zcu.cz, University of West Bohemia
% See paper: TBD doi.

%% Parameters and system simulation
clc
clear variables
close all
format shortG

cd ..
addpath(genpath(pwd)) % Add binaries and data folders

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

    % PF init and param
    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; % initial condition for PF
    ksiPriorsys = ksiPrior; % initial condition for PF
    ksiPriorSTD = ksiPrior; % initial condition for PF
    ksiPriorCor = ksiPrior;


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

    effectiveSampleSize = 0;

    % tic
    % [matSys,matVarSys] = particleFilterEstimationSys(modelChoose, souradniceGNSS, hB, vysky);
    % tocMatSys(mc) = toc;
    % tic
    % [matMult,matVarMult] = particleFilterEstimationMult(modelChoose, souradniceGNSS, hB, vysky);
    % tocMatMult(mc) = toc;

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

        %% Partcle Filter - Multinomial
        tic
        %Measurement Update
        predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); % Prediction density grid through measurement EQ
        predThrMeasEq(1,isnan(predThrMeasEq(1,:))) = inf; % To protect from map extrapolations
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        if effectiveSampleSize < essThrd || k == 1
            w = pdf(V.pdf,pom); % Weights
        else
            w = pdf(V.pdf,pom).*w; % Weights
        end
        w = w/sum(w); % Normalization of weights
        filtMeanPF(:,k) = ksiPrior*w; % mean filtering estimate
        filtVarPF(:,k) = ksiPrior.^2*w - filtMeanPF(:,k).^2; % diagonal of filtering covariance matrix

        covMeas = 0.01 * sqrt(filtVarPF(:, k));  % 1% of the measurement covariance

        effectiveSampleSize = 1/sum(w.^2);

        if effectiveSampleSize < essThrd

            % Resampling
            cumW = cumsum(w);
            randomN = rand(1,noPart);
            % Use this in case the mex file not working
            % for ind = 1:1:noParts
            %     I(ind) = find(cumW >= randomN(ind),1, 'first');
            % end
            I = binarySearch(cumW, randomN, 'first');
            ksi = ksiPrior(:,I);

            % Jitter all particles with multidimensional noise
            jitter = covMeas .* randn(size(ksi));  % Multidimensional jitter
            ksi = ksi + jitter;  % Apply jitter to all particles

        else
            ksi = ksiPrior;
        end

        % Time Update
        ksiPrior = F*ksi + mvnrnd(zeros(1,nx),Q,noPart)' + u(:,k);
        tocPF(k) = toc;


        %% Particle Filter - Systematic
        tic
        %Measurement Update
        predThrMeasEq = hfunct(ksiPriorsys,zeros(nz,1),k); % Prediction density grid through measurement EQ
        predThrMeasEq(1,isnan(predThrMeasEq(1,:))) = inf; % To protect from map extrapolations
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        if effectiveSampleSize < essThrd || k == 1
            w = pdf(V.pdf,pom); % Weights
        else
            w = pdf(V.pdf,pom).*w; % Weights
        end
        w = w/sum(w); % Normalization of weights
        filtMeanPFsys(:,k) = ksiPriorsys*w; % mean filtering estimate
        filtVarPFsys(:,k) = ksiPriorsys.^2*w - filtMeanPFsys(:,k).^2; % diagonal of filtering covariance matrix
        covMeas = 0.01 * sqrt(filtVarPFsys(:, k));  % 1% of the measurement covariance

        effectiveSampleSize = 1/sum(w.^2);

        if effectiveSampleSize < essThrd
            % Resampling
            cumW = cumsum(w);
            randomN = (0:noPart-1)'/noPart + rand/noPart;  % Random starting point + equally spaced
            % Use this in case the mex file not working
            % for ind = 1:1:noParts
            %     I(ind) = find(cumW >= randomN(ind),1, 'first');
            % end
            I = binarySearch(cumW, randomN, 'first');
            ksisys = ksiPriorsys(:,I);

            % Jitter all particles with multidimensional noise
            jitter = covMeas .* randn(size(ksisys));  % Multidimensional jitter
            ksisys = ksisys + jitter;  % Apply jitter to all particles
        else
            ksisys = ksiPriorsys;
        end

        % Time Update
        ksiPriorsys = F*ksisys + mvnrnd(zeros(1,nx),Q,noPart)' + u(:,k);
        tocPFsys(k) = toc;


        %% Particle Filter - STD
        tic
        %Measurement Update
        predThrMeasEq = hfunct(ksiPriorSTD,zeros(nz,1),k); % Prediction density grid through measurement EQ
        predThrMeasEq(1,isnan(predThrMeasEq(1,:))) = inf; % To protect from map extrapolations
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        w = pdf(V.pdf,pom); % Weights
        w = w/sum(w); % Normalization of weights
        filtMeanPFSTD(:,k) = ksiPriorSTD*w; % mean filtering estimate
        filtVarPFSTD(:,k) = ksiPriorSTD.^2*w - filtMeanPFSTD(:,k).^2; % diagonal of filtering covariance matrix

        % Resampling
        cumW = cumsum(w);
        randomN = rand(1,noPart);
        % Use this in case the mex file not working
        % for ind = 1:1:noParts
        %     I(ind) = find(cumW >= randomN(ind),1, 'first');
        % end
        I = binarySearch(cumW, randomN, 'first');
        ksiSTD = ksiPriorSTD(:,I);

        % Time Update
        ksiPriorSTD = F*ksiSTD + mvnrnd(zeros(1,nx),Q,noPart)' + u(:,k);
        tocPFSTD(k) = toc;

        %% Prior-correction particle filter
        tic
        %Measurement Update
        predThrMeasEq = hfunct(ksiPriorCor,zeros(nz,1),k); % Prediction density grid through measurement EQ
        predThrMeasEq(1,isnan(predThrMeasEq(1,:))) = inf; % To protect from map extrapolations
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        if k == 1
            w = pdf(V.pdf,pom); % Weights
        else
            w = pdf(V.pdf,pom).*((exp(sum(-0.5*distVar'*invQ.*distVar',2)))/predDenDenomW); % Weights
        end
        w = w/sum(w); % Normalization of weights
        filtMeanPFCor(:,k) = ksiPriorCor*w; % mean filtering estimate
        filtVarPFCor(:,k) = ksiPriorCor.^2*w - filtMeanPFCor(:,k).^2; % diagonal of filtering covariance matrix

        noPartPrior = floor(0.8*noPart);
        noPartGrid = floor(0.2*noPart);

        alpha = 1.5;

        gridRand = (alpha * sqrt(filtVarPFCor(:,k))) .* rand(nx,noPartGrid);

        % Resampling
        cumW = cumsum(w);
        randomN = rand(1,noPartPrior);
        I = binarySearch(cumW, randomN, 'first');
        ksiCor = ksiPriorCor(:,I);
        ksiCor = [ksiCor gridRand];

        % Time Update
        distVar = mvnrnd(zeros(1,nx),Q,noPart)';
        ksiPriorCor = F*ksiCor + distVar + u(:,k);
        tocPFCor(k) = toc;


        %% Unscented Kalman Filter (basic/recursive/Joseph version)

        tic;
        % UKF - to compare
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


        endtime_act = k;


    end

    % Evaluation for System
    % [rmseMatSys(:, mc), astdMatSys11(mc), astdMatSys22(mc), astdMatSys33(mc), annes_MatSysout(mc)] = ...
    %     evaluateFilter(x, matSys, matVarSys, 'Sys', k); %#ok<*SAGROW>

    % % Evaluation for Mult
    % [rmseMatMult(:, mc), astdMatMult11(mc), astdMatMult22(mc), astdMatMult33(mc), annes_MatMultout(mc)] = ...
    %     evaluateFilter(x, matMult, matVarMult, 'Mult', k);

    % Evaluation for PMF
    [rmsePMF(:, mc), astdPMF11(mc), astdPMF22(mc), astdPMF33(mc), annes_PMFout(mc)] = ...
        evaluateFilter(x, filtMeanPMF, filtVarPMF, 'PMF', k);

    % Evaluation for PFSTD
    [rmsePFSTD(:, mc), astdPF11STD(mc), astdPF22STD(mc), astdPF33STD(mc), annes_PFoutSTD(mc)] = ...
        evaluateFilter(x, filtMeanPFSTD, filtVarPFSTD, 'STD', k);

    % Evaluation for PF
    [rmsePF(:, mc), astdPF11(mc), astdPF22(mc), astdPF33(mc), annes_PFout(mc)] = ...
        evaluateFilter(x, filtMeanPF, filtVarPF, 'PF', k);

    % Evaluation for PFCor
    [rmsePFCor(:, mc), astdPF11Cor(mc), astdPF22Cor(mc), astdPF33Cor(mc), annes_PFoutCor(mc)] = ...
        evaluateFilter(x, filtMeanPFCor, filtVarPFCor, 'Cor', k);

    % Evaluation for PFsys
    [rmsePFsys(:, mc), astdPF11sys(mc), astdPF22sys(mc), astdPF33sys(mc), annes_PFoutsys(mc)] = ...
        evaluateFilter(x, filtMeanPFsys, filtVarPFsys, 'sys', k);

    % Evaluation for UKF
    [rmseUKF(:, mc), astdUKF11(mc), astdUKF22(mc), astdUKF33(mc), annes_UKFout(mc)] = ...
        evaluateFilter(x, xf_ukf, Pf_ukf, 'UKF', k);

    tocPMFavg(:,mc) = mean(tocPMF);
    tocPFavgCor(mc) = mean(tocPFCor);
    tocPMFavg(mc) =  mean(tocPMF);
    tocPFavg(mc) =  mean(tocPF);
    tocPFavgSTD(mc) =  mean(tocPFSTD);
    tocPFavgsys(mc) =  mean(tocPFsys);
    tocUKFavg(mc) =  mean(tocUKF);

end

% % Normalize outputs for System and Mult (commented out)
% annes_MatSysout = annes_MatSysout / (nx * k);
% annes_MatMultout = annes_MatMultout / (nx * k);
annes_PMFout = annes_PMFout / (nx * k);
annes_PFout = annes_PFout / (nx * k);
annes_PFoutsys = annes_PFoutsys / (nx * k);
annes_PFoutSTD = annes_PFoutSTD / (nx * k);
annes_UKFout = annes_UKFout / (nx * k);

rmsePMFout = mean(rmsePMF, 2);
rmsePFout = mean(rmsePF, 2);
rmsePFsysout = mean(rmsePFsys, 2);
rmsePFSTDout = mean(rmsePFSTD, 2);
rmseUKFout = mean(rmseUKF, 2);

tocPMFavgOut = mean(tocPMFavg, 2);
tocPFavgOut = mean(tocPFavg, 2);
tocPFavgSTDOut = mean(tocPFavgSTD, 2);
tocPFavgsysOut = mean(tocPFavgsys, 2);
tocUKFavgOut = mean(tocUKFavg, 2);
% tocMatavgSysOut = mean(tocMatSys/(endTime-1), 2); % Commented out
% tocMatavgMultOut = mean(tocMatMult/(endTime-1), 2); % Commented out
tocPFavgCorOut = mean(tocPFavgCor, 2);

% % Normalize outputs for System and Mult (commented out)
% annes_MatSysout = annes_MatSysout / (nx * k);
% annes_MatMultout = annes_MatMultout / (nx * k);
annes_PFCorout = annes_PFoutCor / (nx * k);

% % Mean RMSE calculations for System and Mult (commented out)
% rmseMatSysout = mean(rmseMatSys, 2);
% rmseMatMultout = mean(rmseMatMult, 2);
rmsePFCorout = mean(rmsePFCor, 2);  % Add PF Cor RMSE calculation

% % Mean standard deviations for System and Mult (commented out)
% astdMatSys = [mean(astdMatSys11), mean(astdMatSys22)];
% astdMatMult = [mean(astdMatMult11), mean(astdMatMult22)];
astdPFCor = [mean(astdPF11Cor), mean(astdPF22Cor)];  % Add PF Cor standard deviations

meanTIME = [mean(tocPFavgCorOut), mean(tocPMFavgOut), mean(tocPFavgOut), ...
    mean(tocPFavgsysOut), mean(tocPFavgSTDOut), mean(tocUKFavgOut)]';

% Prepare other required mean values for table
meanRMSEs = [mean(rmsePFCorout(1:2)), mean(rmsePMFout(1:2)), mean(rmsePFout(1:2)), ...
    mean(rmsePFsysout(1:2)), mean(rmsePFSTDout(1:2)), mean(rmseUKFout(1:2))]';

meanASTDs = [mean(astdPFCor), mean([astdPMF11 astdPMF22]), mean([astdPF11 astdPF22]), ...
    mean([astdPF11sys astdPF22sys]), mean([astdPF11STD astdPF22STD]), ...
    mean([astdUKF11 astdUKF22])]';

meanANNES = [mean(annes_PFCorout), mean(annes_PMFout), mean(annes_PFout), ...
    mean(annes_PFoutsys), mean(annes_PFoutSTD), mean(annes_UKFout)]';

% Create table T2 including all filtering methods
T2 = table(...
    meanRMSEs, ...
    meanASTDs, ...
    meanANNES, ...
    meanTIME, ...
    'VariableNames', {'RMSE POS','ASTD POS','ANNES','TIME'}, ...
    'RowName', {'PF Cor', 'PMF','PF multinomial, ess, jitter','PF systematic, ess, jitter','PF std multi','UKF'} ...
    );