%% Lagrangian Grid-Based Filter (LGbF) and others 
%   for real TAN dataset
%   author: pesslovany@gmail.com

clc
clear variables
close all
format shortG

addpath(genpath(pwd)); % Add all files and subfolders in the current directory

% Select Filters to Run
runFlags.LGF_Standard = true;     % classic LGbF
runFlags.LGbF_Spectral = false;     % spectral LGbF
runFlags.RBPF = false;           % Rao-Blackewellized Particle Filter
runFlags.PF = false;           % Bootstrap Particle Filter
runFlags.UKF = false;           % Unscented Kalman Filter
runFlags.EnGMF = true;          % Ensemble Gaussian Mixture Filter

% Parameters and system simulation
modelChoose = 4; % choose model: 3D = 3, 4D = 4

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

    [predGrid, predGridDelta, gridDimOld, xOld, Ppold, ~] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    [predGrid2, predGridDelta2, gridDimOld2, xOld2, Ppold2, ~] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    pom = (predGrid - meanX0);
    denominator = sqrt((2*pi)^nx * det(varX0));
    predDensityProb  = exp(sum(-0.5 * pom'/(varX0) .* pom',2)) / denominator;
    predDensityProb2 = predDensityProb;
    predPdf3         = predDensityProb;

    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; % initial condition for PF

    ensemblePosterior = mvnrnd(meanX0,varX0,noEns)'; % initial condition for EnGMF

    % UKF init
    xp = meanX0;
    Pp = varX0;

    predDenDenomW = sqrt((2*pi)^nx * det(Q));
    halfGrid = ceil(N/2);
    GridDelta2 = [];
    GridDelta = [];

    predW = 1/noPartRbpf*ones(noPartRbpf,1);
    predXn = meanX0(xnid) + chol(varX0(xnid,xnid),'lower')*randn(nxn, noPartRbpf); % assumes that the initial nonlinear and linear states are independent
    predXl = repmat(meanX0(xlid), 1, noPartRbpf);
    predPl = repmat(varX0(xlid,xlid), 1, 1, noPartRbpf);

    measMeanEnGMF(:,1) = meanX0;
    covEnGMF(:,:,1) = varX0;

    for k = 1:endTime-1
        disp(['Step:', num2str(k), '/', num2str(endTime-1)])

        %% Ensemble- Gaussian Mixture Filter
        if runFlags.EnGMF
        tic

        wbark = ones(noEns,1);
        wbark = wbark/(sum(wbark));
        Xbark = F*ensemblePosterior;
        Ps = diracToGaussMix(Xbark,wbark,Q);

        [XkGSF,PkGSF,wkGSF] = ukfGSMupdate(Xbark,wbark,nx,Ps,z,k,hfunct,R);
        
        measMeanEnGMF(:,k+1) = XkGSF * wkGSF;
        covEnGMF(:,:,k+1) = sum(PkGSF .* reshape(wkGSF,1,1,[]), 3);
        nuxk  = XkGSF - measMeanEnGMF(:,k+1);
        covEnGMF(:,:,k+1) = covEnGMF(:,:,k+1) + (nuxk .* repmat(wkGSF',nx,1)) * nuxk';
        covEnGMF(:,:,k+1) = (covEnGMF(:,:,k+1) + covEnGMF(:,:,k+1).')/2;

        resampleInd = systematicResampling(wkGSF,noEns);

        [ensemblePosterior] = realisationOfGauss(resampleInd,XkGSF,PkGSF);

        tocEnGMF(k) = toc; 
        end

        %% Standard GBF
        if runFlags.LGF_Standard
            tic
            [measPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predDensityProb,predGridDelta(:,k),hfunct);
            [measMeanGBF(:,k), measVarGBF(:,:,k)] = momentGbF(predGrid, measPdf, predGridDelta(:,k));
            measGrid = predGrid;

            [measPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(measPdf, ...
                gridDimOld,xOld,Ppold,Q,sFactor,nx,Npa,measMeanGBF(:,k),measVarGBF(:,:,k),F);
            xOld = F * measMeanGBF(:,k) + u(:,k);
            Ppold = F * eigVect;
            [predDensityProb,predGrid,predGridDelta] = gbfTimeUpdateFFT(F,measPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));

            tocPMF(k) = toc;
        end

        %% Spectral GBF
        if runFlags.LGbF_Spectral
            tic
            [measPdf2] = gbfMeasUpdate(predGrid2,nz,k,z(:,k),V,predDensityProb2,predGridDelta2(:,k),hfunct);
            [measMeanSGBF(:,k), measVarSGBF(:,:,k)] = momentGbF(predGrid2, measPdf2, predGridDelta2(:,k));
            measGrid2 = predGrid2;

            gridCenter = xOld2;
            gridRotation = Ppold2;
            gridDim = gridDimOld2;

            [measPdf2, GridDelta2, xOld2, predGridAdvect, gridDimOld2, ...
                Ppold2, gridBound, measGridNew2, predVarPom] = ...
                measPdfInterpSpect(measPdf2, GridDelta2, measMeanSGBF, measVarSGBF, ...
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
            [measMeanRBPF(:,k), measVarRBPF(:,:,k)] = momentRBPF(measXn, measXl, measPl, measW, nxn, nxl, xnid, xlid);

            [predXn, predXl, predPl, predW] = rbpfTimeUpdate(measXn, measXl, measPl, measW, F, Q, u(:,k), nxn, nxl, xnid, xlid);
            tocRBPF(k) = toc;
        end

        %% Bootstrap PF
        if runFlags.PF
            tic
            %Measurement Update
            predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); % Prediction density grid through measurement EQ
            predThrMeasEq(1,isnan(predThrMeasEq(1,:))) = inf; % To protect from map extrapolations
            pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
            w = pdf(V.pdf,pom); % Weights
            w = w/sum(w); % Normalization of weights
            measMeanPF(:,k) = ksiPrior*w; % mean filtering estimate
            measVarPF(:,:,k) = diag(ksiPrior.^2*w - measMeanPF(:,k).^2); % diagonal of filtering covariance matrix

            [resampleInd] = systematicResampling(w,noPart);
            ksi = ksiPrior(:, resampleInd);

            % Time Update
            ksiPrior = F*ksi + mvnrnd(zeros(1,nx),Q,noPart)' + u(:,k);
            tocPF(k) = toc;
        end

        %% Unscented Kalman filter
        if runFlags.UKF
            tic
            % Meas update
            [chi, wm, wc] = msp(xp, Pp, kappa);
            dzetap = hfunct(chi, zeros(nz,1), k);
            zp = dzetap * wm' + meanV;
            dz = dzetap - zp;
            Pzp = dz * diag(wc) * dz' + R;
            dx = chi - xp;
            Pxzp = dx * diag(wc) * dz';
            K = Pxzp / Pzp;
            e = z(:,k) - zp;
            measMeanUKF(:,k) = xp + K*e;
            measVarUKF(:,:,k) = Pp - K*Pzp*K';
            % Time update
            xp = F*measMeanUKF(:,k) + u(:,k);
            Pp = F*measVarUKF(:,:,k)*F' + Q;
            tocUKF(k) = toc;
        end

    end

    % Evaluation
    if runFlags.LGF_Standard
        [rmsePMF(:,mc), astdPMF(:,mc), annesPMF(mc)] = calcRes(x, measMeanGBF, measVarGBF, k); %#ok<*SAGROW>
        tocPMFavg(mc) = mean(tocPMF);
    end
    if runFlags.LGbF_Spectral
        [rmsePMF2(:,mc), astdPMF2(:,mc), annesPMF2(mc)] = calcRes(x, measMeanSGBF, measVarSGBF, k);
        tocPMF2avg(mc) = mean(tocPMF2);
    end
    if runFlags.RBPF
        [rmseRBPF(:,mc), astdRBPF(:,mc), annesRBPF(mc)] = calcRes(x, measMeanRBPF, measVarRBPF, k);
        tocRBPFavg(mc) = mean(tocRBPF);
    end
    if runFlags.PF
        [rmsePF(:,mc), astdPF(:,mc), annesPF(mc)] = calcRes(x, measMeanPF, measVarPF, k);
    end
    if runFlags.UKF
        [rmseUKF(:,mc), astdUKF(:,mc), annesUKF(mc)] = calcRes(x, measMeanUKF, measVarUKF, k);
    end
    if runFlags.EnGMF
        [rmseEnGMF(:,mc), astdEnGMF(:,mc), annesEnGMF(mc)] = calcRes(x, measMeanEnGMF, covEnGMF, k);
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
if runFlags.PF
    annes_PFoutPom = sum(annesPF) / (nx * mc * (k + 1));
    rmsePFout      = mean(rmsePF, 2);
    astdPFout      = mean(astdPF, 2);
    tocPFavgout    = mean(tocPF);
end
if runFlags.UKF
    annes_UKFoutPom = sum(annesUKF) / (nx * mc * (k + 1));
    rmseUKFout      = mean(rmseUKF, 2);
    astdUKFout      = mean(astdUKF, 2);
    tocUKFavgout    = mean(tocUKF);
end
if runFlags.EnGMF
    annes_EnGMFoutPom = sum(annesEnGMF) / (nx * mc * (k + 1));
    rmseEnGMFout      = mean(rmseEnGMF, 2);
    astdEnGMFout      = mean(astdEnGMF, 2);
    tocEnGMFavgout    = mean(tocEnGMF);
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
if runFlags.PF
    rowNames{end+1} = 'PF bootstrap';
    vals{end+1} = [rmsePFout(1), rmsePFout(2), rmsePFout(3), ...
        astdPFout(1), astdPFout(2), astdPFout(3), ...
        annes_PFoutPom, tocPFavgout];
end
if runFlags.UKF
    rowNames{end+1} = 'UKF';
    vals{end+1} = [rmseUKFout(1), rmseUKFout(2), rmseUKFout(3), ...
        astdUKFout(1), astdUKFout(2), astdUKFout(3), ...
        annes_UKFoutPom, tocUKFavgout];
end
if runFlags.EnGMF
    rowNames{end+1} = 'En-GMF';
    vals{end+1} = [rmseEnGMFout(1), rmseEnGMFout(2), rmseEnGMFout(3), ...
        astdEnGMFout(1), astdEnGMFout(2), astdEnGMFout(3), ...
        annes_EnGMFoutPom, tocEnGMFavgout];
end

T2 = cell2mat(vals');
T2 = array2table(T2, ...
    'VariableNames', {'RMSE_x1','RMSE_x2','RMSE_x3','ASTD_1','ASTD_2','ASTD_3','ANNES','TIME'}, ...
    'RowNames', rowNames);

disp(T2)

% Trajectory Plots
plot(x(1,:),x(2,:)); hold on;
if runFlags.LGF_Standard, plot(measMeanGBF(1,:),measMeanGBF(2,:)); end
if runFlags.LGbF_Spectral, plot(measMeanSGBF(1,:),measMeanSGBF(2,:)); end
if runFlags.RBPF, plot(measMeanRBPF(1,:),measMeanRBPF(2,:)); end
if runFlags.PF, plot(measMeanPF(1,:),measMeanPF(2,:)); end
if runFlags.UKF, plot(measMeanUKF(1,:),measMeanUKF(2,:)); end
if runFlags.EnGMF, plot(measMeanEnGMF(1,:),measMeanEnGMF(2,:)); end


