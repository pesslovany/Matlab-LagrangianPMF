%%
%   Efficient Point-Mass Filter for 4D nearly constant velocity using CPD tensor decomposition.
%   author: pesslovany@gmail.com
%   Papers reference will be added after publication
%% Parameters and system simulation

warning('off','MATLAB:nearlySingularMatrix'); % CAREFULL
addpath(genpath(pwd)); % Add all files and subfolders in the current directory


% clc
clear
% rng(1);
close all
format shortG
load('data.mat') % map of terrain

% souradniceGNSS = souradniceGNSS(:,1:50);
% hB = hB(1:50);
vysky = griddedInterpolant(souradniceX',souradniceY',souradniceZ',"linear","nearest");
timeSteps = 1:1:length(souradniceGNSS);

% Max rank of CPD
Rdef = 20;
% MC simulations
MC = 10;

for mc = 1:1:MC

    model = initModelCPD(souradniceGNSS, hB, vysky);
    fields = fieldnames(model);
    for i = 1:numel(fields)
        eval([fields{i} ' = model.' fields{i} ';']);
    end
    clear model

    Q = diag(diag(Q));
    

    % PMF init and param
    % Initial grid
    [predGrid, predGridDelta, gridDimOld, xOld, Ppold] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    [predGrid2, predGridDelta2, gridDimOld2, xOld2, Ppold2] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial PMD
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)
    % Auxiliary variables
    predDenDenomW = sqrt((2*pi)^nx*det(Q)); % Denominator for convolution in predictive step
    halfGrid = ceil(N/2); % Middle row of the TPM matrix index

    % Intial CPD PMD
    denominator1D = sqrt((2*pi)*diag(varX0));
    gridCPD = arrayfun(@(m, v, nP) linspace(m - sFactorCPD*sqrt(v), m + sFactorCPD*sqrt(v), nP), ...
        meanX0, diag(varX0), NpaCPD', 'UniformOutput', false);
    predGridDeltaCPD = cellfun(@(x) x(2)-x(1), gridCPD);
    for indDim = 1:nx
        pom = gridCPD{indDim} - meanX0(indDim);
        gridDimCPD{indDim} = ((exp(-0.5*pom/(varX0(indDim,indDim)).*pom))/denominator1D(indDim))';
    end
    predDensityProbCPD = ktensor(gridDimCPD);


    for k = 1:1:endTime-1
        disp(['Step:', num2str(k),'/',num2str(endTime-1)])
        %% Grid-based Filter
        tic

        % Measurement update
        [measPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predDensityProb,predGridDelta(:,k),hfunct); % Measurement update
        % measPdf(measPdf<1e-11) = 0;
        % Filtering mean and var
        measMean2(:,k) = predGrid*measPdf*prod(predGridDelta(:,k)); % Measurement update mean
        chip_ = (predGrid-measMean2(:,k));
        chip_w = chip_.*repmat(measPdf',nx,1);
        measVar2(:,:,k) = chip_w*chip_' * prod(predGridDelta(:,k)); % Measurement update variance
        % measVaricka = diag(measVar2(:,:,k))
        measGrid = predGrid;
        [measPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(measPdf,...
            gridDimOld,xOld,Ppold,Q,sFactor,nx,Npa,measMean2(:,k),measVar2(:,:,k),F);
        xOld = F*measMean2(:,k) + u(:,k);
        Ppold = F*(eigVect);
        % Time Update
        [predDensityProb,predGrid,predGridDelta] = gbfTimeUpdateFFT(F,measPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));
        % predDensityProb(predDensityProb<0)=0;
        tocPMF(k) = toc; % Time evaluation

        %% RES
        % meas = predGrid * predDensityProb * prod(predGridDelta(:,k+1));
        % chip_ = predGrid - meas;
        % chip_w = chip_ .* repmat(predDensityProb', nx, 1);
        % measV = chip_w * chip_' * prod(predGridDelta(:,k+1));
        % 
        % fprintf('\n%-12s | %-12s | %-12s | %-12s\n', ...
        %     'meas', 'diag(meas)', 'pred', 'diag(pred)');
        % fprintf('%s\n', repmat('-', 1, 54));
        % for i = 1:nx
        %     fprintf('%12.4f | %12.4f | %12.4f | %12.4f\n', ...
        %         measMean2(i,k), ...
        %         measVar2(i,i,k), ...
        %         meas(i), ...
        %         measV(i,i));
        % end

        %% Grid-based Filter - cont - CPD

        tic
        [measPdfCPD] = gbfMeasUpdateCPD( ...
            gridCPD, z(:,k), Vsplit, ...
            predDensityProbCPD, predGridDeltaCPD(:,k), ...
            hFunctSplit, mapping, NpaCPD, Rdef);
        [measMeanCPD(:,k), measVarCPD(:,k)] = ...
            momentsCPD(measPdfCPD, gridCPD, predGridDeltaCPD(:,k));

        % Advection 
        [predDensityCPD, gridCPD, predGridDeltaCPD] = gbfAdvectionCPD(measPdfCPD, measMeanCPD(:,k), measVarCPD(:,k), ...
                                                    gridCPD, F, Q, NpaCPD, sFactorCPD, predGridDeltaCPD, Rdef, k);

        % [measPredAdv, measVPredAdv] = momentsCPD(predDensityCPD, gridCPD, predGridDeltaCPD(:,k+1));

        % Spectral diffusion
        [predDensityProbCPD] = gbfDiffusionCPD(predDensityCPD, Q, predGridDeltaCPD, nx, k, Rdef);
        % Discrete diffusion
        % [predDensityProbCPD] = gbfDiffusionCPDSpect(predDensityCPD, Q, predGridDeltaCPD, nx, k, Rdef);

        tocPMFCPD(k) = toc; % Time evaluation

        % plot(reshape(double(full(measPdfCPD)),1,[]))
        
        %% Res
        % [measPred, measVPred] = momentsCPD(predDensityProbCPD, gridCPD, predGridDeltaCPD(:,k+1));
        % %--- trimmed 4-column table ---
        % fprintf('\n%-12s | %-12s | %-12s | %-12s\n', ...
        %     'measCPD', 'diag(measCPD)', 'predCPD', 'diag(predCPD)');
        % fprintf('%s\n', repmat('-', 1, 51));   % adjust width to match header
        % for i = 1:nx
        %     fprintf('%12.4f | %12.4f | %12.4f | %12.4f\n', ...
        %         measMeanCPD(i,k), ...
        %         measVarCPD(i,k), ...
        %         measPred(i), ...
        %         measVPred(i));
        % end
        % 
        % trueAdvCPD = F*diag(measVarCPD(:,k))*F';
        % 
        % fprintf('\n%-12s | %-12s | %-12s\n', ...
        %     'trueAdvCPD', 'diag(advCPD)', 'advCPD');
        % % total width = 12 + 3 + 12 + 3 + 12 = 42
        % fprintf('%s\n', repmat('-', 1, 42));
        % for i = 1:nx
        %     fprintf('%12.4f | %12.4f | %12.4f\n', ...
        %         trueAdvCPD(i,i), measVPredAdv(i), ...
        %         measPredAdv(i));
        % end

    end

    % k=k+1;
    % Evaluation
    rmsePMF(:,mc) = sqrt(mean((x(:,1:k-1)-measMean2(:,1:k-1)).^2,2)); %#ok<*SAGROW>
    astdPMF11(mc) = sqrt(mean(measVar2(1,1,1:k-1)));
    astdPMF22(mc) = sqrt(mean(measVar2(2,2,1:k-1)));
    astdPMF33(mc) = sqrt(mean(measVar2(3,3,1:k-1)));
    astdPMF44(mc) = sqrt(mean(measVar2(4,4,1:k-1)));
    annes_PMF = 0;
    for indAn = 1:1:k-1
        annes_PMF = annes_PMF + ((x(:,indAn)-measMean2(:,indAn)).*(1./(diag(measVar2(:,:,indAn)))))'*(x(:,indAn)-measMean2(:,indAn));
    end
    annes_PMFout(mc) = annes_PMF;


    % Evaluation
    rmsePMFCPD(:,mc) = sqrt(mean((x(:,1:k-1)-measMeanCPD(:,1:k-1)).^2,2));
    astdPMF11CPD(mc) = sqrt(mean(measVarCPD(1,1:k-1)));
    astdPMF22CPD(mc) = sqrt(mean(measVarCPD(2,1:k-1)));
    astdPMF33CPD(mc) = sqrt(mean(measVarCPD(3,1:k-1)));
    astdPMF44CPD(mc) = sqrt(mean(measVarCPD(4,1:k-1)));
    annes_PMFCPD = 0;
    for indAn = 1:1:k-1
        annes_PMFCPD = annes_PMFCPD + ((x(:,indAn)-measMeanCPD(:,indAn)).*(1./((measVarCPD(:,indAn)))))'*(x(:,indAn)-measMeanCPD(:,indAn));
    end
    annes_PMFoutCPD(mc) = annes_PMFCPD;


    tocPMFavgCPD(mc) = mean(tocPMFCPD);
    tocPMFavg(mc) = mean(tocPMF);
end

annes_PMFout = sum(annes_PMFout)/(nx*mc*(k+1));
annes_PMFoutCPD = sum(annes_PMFCPD)/(nx*mc*(k+1));

rmsePMFout = mean(rmsePMF,2);
rmsePMFoutCPD = mean(rmsePMFCPD,2);
tocPMFavgOut = mean(tocPMFavg,2);

T2 = table([ rmsePMFout(1)  rmsePMFoutCPD(1)]',...
    [ rmsePMFout(2)  rmsePMFoutCPD(2)]',...
    [ rmsePMFout(3)  rmsePMFoutCPD(3)]',...
    [ rmsePMFout(4) rmsePMFoutCPD(4)]',...
    [ mean(astdPMF11)  real(mean(astdPMF11CPD))]',...
    [ mean(astdPMF22) real(mean(astdPMF22CPD))]',...
    [ mean(astdPMF33) real(mean(astdPMF33CPD))]',...
    [ mean(astdPMF44) real(mean(astdPMF44CPD))]',...
    [ mean(annes_PMFout) real(mean(annes_PMFoutCPD))]',...
    [ mean(tocPMFavgOut) mean(tocPMFavgCPD)]',...
    'VariableNames',{'RMSE x1','RMSE x2','RMSE x3','RMSE x4','ASTD 1','ASTD 2','ASTD 3','ASTD 4','ANNES','TIME'},...
    'RowName',{'PMF','CPD PMF'}) %#ok<NOPTS>


plot(x(1,:),x(2,:))
hold on
plot(measMean2(1,:),measMean2(2,:))
plot(measMeanCPD(1,:),measMeanCPD(2,:))