%% Lagrangian Grid-Based Filter (LGbF) Setup
%   FFT and spectral differentiation based filters for real TAN dataset
%   author: pesslovany@gmail.com

clc
clear variables
close all
format shortG

addpath(genpath(pwd)); % Add all files and subfolders in the current directory

% Select Filters to Run
runFlags.LGF_Standard = true;     % classic LGbF
runFlags.LGbF_Spectral = true;     % spectral LGbF
runFlags.RBPF = true;           % Rao-Blackewellized Particle Filter

% Parameters and system simulation
modelChoose = 3; % choose model: 3D = 3, 4D = 4

load('data.mat') % map of terrain
vysky = griddedInterpolant(souradniceX',souradniceY',souradniceZ',"linear","none");

timeSteps = 1:1:length(souradniceGNSS);
endTime = length(timeSteps);
MC = 1;

for mc = 1:1:MC
    model = initModel(modelChoose, souradniceGNSS, hB, vysky);
    fields = fieldnames(model);
    for i = 1:numel(fields)
        eval([fields{i} ' = model.' fields{i} ';']);
    end
    clear model

    sFactor = 6;

    [predGrid, predGridDelta, gridDimOld, xOld, Ppold, ~] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    [predGrid2, predGridDelta2, gridDimOld2, xOld2, Ppold2, ~] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    pom = (predGrid - meanX0);
    denominator = sqrt((2*pi)^nx * det(varX0));
    predDensityProb  = exp(sum(-0.5 * pom'/(varX0) .* pom',2)) / denominator;
    predDensityProb2 = predDensityProb;
    predPdf3         = predDensityProb;

    predDenDenomW = sqrt((2*pi)^nx * det(Q));
    halfGrid = ceil(N/2);
    GridDelta2 = [];
    GridDelta = [];

    predW = 1/noPartRbpf*ones(noPartRbpf,1);
    predXn = meanX0(xnid) + chol(varX0(xnid,xnid),'lower')*randn(nxn, noPartRbpf); % assumes that the initial nonlinear and linear states are independent
    predXl = repmat(meanX0(xlid), 1, noPartRbpf);
    predPl = repmat(varX0(xlid,xlid), 1, 1, noPartRbpf);

    for k = 1:endTime-1
        disp(['Step:', num2str(k), '/', num2str(endTime-1)])

        %% Standard GBF
        if runFlags.LGF_Standard
            tic
            [measPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predDensityProb,predGridDelta(:,k),hfunct);
            [measMean2(:,k), measVar2(:,:,k)] = momentGbF(predGrid, measPdf, predGridDelta(:,k));
            measGrid = predGrid;

            [measPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(measPdf, ...
                gridDimOld,xOld,Ppold,Q,sFactor,nx,Npa,measMean2(:,k),measVar2(:,:,k),F);
            xOld = F * measMean2(:,k) + u(:,k);
            Ppold = F * eigVect;
            [predDensityProb,predGrid,predGridDelta] = gbfTimeUpdateFFT(F,measPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));

            tocPMF(k) = toc;
        end

        %% Spectral GBF
        if runFlags.LGbF_Spectral
            tic
            [measPdf2] = gbfMeasUpdate(predGrid2,nz,k,z(:,k),V,predDensityProb2,predGridDelta2(:,k),hfunct);
            [measMean3(:,k), measVar3(:,:,k)] = momentGbF(predGrid2, measPdf2, predGridDelta2(:,k));
            measGrid2 = predGrid2;

            gridCenter = xOld2;
            gridRotation = Ppold2;
            gridDim = gridDimOld2;

            [measPdf2, GridDelta2, xOld2, predGridAdvect, gridDimOld2, ...
                Ppold2, gridBound, measGridNew2, predVarPom] = ...
                measPdfInterpSpect(measPdf2, GridDelta2, measMean3, measVar3, ...
                k, u, F, Q, sFactor, nx, Npa, ...
                gridDim, gridRotation, gridCenter);

            rotQ = Ppold2' * Q * Ppold2;
            [predDensityProb2,predGrid2,predGridDelta2] = gbfTimeUpdateSpect(F,measPdf2,measGridNew2,GridDelta2,k,Npa,rotQ,u(:,k),gridBound);
            if abs(min(predDensityProb2)) > max(predDensityProb2)/100
                [predDensityProb2,predGrid2,predGridDelta2] = gbfTimeUpdateFFT(F,measPdf2,measGridNew2,GridDelta2,k,Npa,invQ,predDenDenomW,nx,u(:,k));
            end
            predDensityProb2(predDensityProb2 < 0) = 0;
            predDensityProb2 = predDensityProb2 / (sum(predDensityProb2)*prod(GridDelta2(:,k+1)))';
            tocPMF2(k) = toc;
        end
        %% RBPF
        if runFlags.RBPF
            tic
            [measXn, measXl, measPl, measW] = rbpfMeasUpdate(predXn, predXl, predPl, predW, nx, nxl, nz, xnid, xlid, k, z(:,k), R, hfunct, hJacobLfunct, essThrdRbpf);
            [measMean4(:,k), measVar4(:,:,k)] = momentRBPF(measXn, measXl, measPl, measW, nxn, nxl, xnid, xlid);

            [predXn, predXl, predPl, predW] = rbpfTimeUpdate(measXn, measXl, measPl, measW, F, Q, u(:,k), nxn, nxl, xnid, xlid);
            tocRBPF(k) = toc;
        end

    end

    % Evaluation
    if runFlags.LGF_Standard
        [rmsePMF(:,mc), astdPMF(:,mc), annesPMF(mc)] = calcRes(x, measMean2, measVar2, k); %#ok<*SAGROW>
        tocPMFavg(mc) = mean(tocPMF);
    end
    if runFlags.LGbF_Spectral
        [rmsePMF2(:,mc), astdPMF2(:,mc), annesPMF2(mc)] = calcRes(x, measMean3, measVar3, k);
        tocPMF2avg(mc) = mean(tocPMF2);
    end
    if runFlags.RBPF
        [rmseRBPF(:,mc), astdRBPF(:,mc), annesRBPF(mc)] = calcRes(x, measMean4, measVar4, k);
        tocRBPFavg(mc) = mean(tocRBPF);
    end
end

% Post-Processed Statistics
if runFlags.LGF_Standard
    annes_PMFoutPom  = sum(annesPMF)  / (nx * mc * (k + 1));
    rmsePMFout       = mean(rmsePMF , 2);
    astdPMFout       = mean(astdPMF , 2);
    tocPMFavgOut     = mean(tocPMFavg);
end
if runFlags.LGbF_Spectral
    annes_PMFout2Pom = sum(annesPMF2) / (nx * mc * (k + 1));
    rmsePMFout2      = mean(rmsePMF2, 2);
    astdPMFout2      = mean(astdPMF2, 2);
    tocPMFavg2out    = mean(tocPMF2avg);
end
if runFlags.RBPF
    annes_RBPFoutPom = sum(annesRBPF) / (nx * mc * (k + 1));
    rmseRBPFout      = mean(rmseRBPF, 2);
    astdRBPFout      = mean(astdRBPF, 2);
    tocRBPFavgout    = mean(tocRBPFavg);
end

% Final Results Table
rowNames = {};
vals = {};
if runFlags.LGF_Standard
    rowNames{end+1} = 'LGbF';
    vals{end+1} = [rmsePMFout(1), rmsePMFout(2), rmsePMFout(3), ...
        astdPMFout(1), astdPMFout(2), astdPMFout(3), ...
        annes_PMFoutPom, tocPMFavgOut];
end
if runFlags.LGbF_Spectral
    rowNames{end+1} = 'Spect LGbF';
    vals{end+1} = [rmsePMFout2(1), rmsePMFout2(2), rmsePMFout2(3), ...
        astdPMFout2(1), astdPMFout2(2), astdPMFout2(3), ...
        annes_PMFout2Pom, tocPMFavg2out];
end
if runFlags.RBPF
    rowNames{end+1} = 'RBPF';
    vals{end+1} = [rmseRBPFout(1), rmseRBPFout(2), rmseRBPFout(3), ...
        astdRBPFout(1), astdRBPFout(2), astdRBPFout(3), ...
        annes_RBPFoutPom, tocRBPFavgout];
end

T2 = cell2mat(vals');
T2 = array2table(T2, ...
    'VariableNames', {'RMSE_x1','RMSE_x2','RMSE_x3','ASTD_1','ASTD_2','ASTD_3','ANNES','TIME'}, ...
    'RowNames', rowNames);

disp(T2)

% Trajectory Plots
if runFlags.LGF_Standard, plot(x(1,:),x(2,:)); hold on; plot(measMean2(1,:),measMean2(2,:)); end
if runFlags.LGbF_Spectral, plot(measMean3(1,:),measMean3(2,:)); end
if runFlags.RBPF, plot(measMean4(1,:),measMean4(2,:)); end
