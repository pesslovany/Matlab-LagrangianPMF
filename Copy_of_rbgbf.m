%% Lagrangian Grid-Based Filter (LGbF) Setup
%   FFT and spectral differentiation based filters for real TAN dataset
%   author: pesslovany@gmail.com

clc
clear variables
close all
format shortG
%rng(1)
addpath(genpath(pwd)); % Add all files and subfolders in the current directory

% Select Filters to Run
runFlags.LGF_Standard = true;     % classic LGbF
runFlags.LGbF_Spectral = true;     % spectral LGbF

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
        predDenDenomW2 = sqrt((2*pi)^2 * det(Q(1:2,1:2)));

    halfGrid = ceil(N/2);
    GridDelta2 = [];
    GridDelta = [];

    % RBPMF matrices initialisation (MODEL 3 ONLY considered)
    % --
    nn = 2;
    nl = 1;
    F11 = F(1:nn,1:nn);
    F12 = F(1:nn,nn+1:nn+nl);
    F21 = F(nn+1:nn+nl,1:nn);
    F22 = F(nn+1:nn+nl,nn+1:nn+nl);
    Q11 = Q(1:nn,1:nn);
    Q12 = Q(1:nn,nn+1:nn+nl);
    Q21 = Q(nn+1:nn+nl,1:nn);
    Q22 = Q(nn+1:nn+nl,nn+1:nn+nl);
    meanX01 = meanX0(1:nn);
    meanX02 = meanX0(nn+1:nn+nl);
    Pinit11 = varX0(1:nn,1:nn);
    Pinit22 = varX0(nn+1:nn+nl,nn+1:nn+nl);
    Pinit = blkdiag(Pinit11, Pinit22);
    PMG_max_r = Npa(1);
    PMG_max_c = Npa(2);
    % -- measurement function def
    h1 =  @(xn) vysky(xn(1,:),xn(2,:));
    H2 = [1];
    invR = inv(R);
    % rNorm = 1/(2*pi*sqrt(det(Qaux)));
    mvnpdf1 = @(x, mu, P, invP, n) 1/sqrt((2*pi)^n*det(P)) * exp(-0.5*(x-mu)'*invP*(x-mu));
    mvnpdf1diag2D = @(x, mu, P, invP) 1/(2*pi*sqrt(P(1,1)*P(2,2))) * exp(-0.5*( (x(1,:)-mu(1)).^2*invP(1,1) + (x(2,:)-mu(2)).^2*invP(2,2)));
    % -- state estimate initialisation
    Nrb = PMG_max_r * PMG_max_c;
    idxVec = (1:Nrb)';
    initSigmaCoverage = sFactor;
    Cov_PMG0_m = Pinit(1:2,1:2);
    % std_PMG0_m11 = sqrt(Cov_PMG0_m(1,1));
    % std_PMG0_m22 = sqrt(Cov_PMG0_m(2,2));
    % % - grid computation
    % PMG_Evec = linspace(-initSigmaCoverage*std_PMG0_m11,initSigmaCoverage*std_PMG0_m11,PMG_max_c)+meanX01(1);
    % PMG_deltaE = PMG_Evec(2)-PMG_Evec(1);
    % PMG_Nvec = linspace(-initSigmaCoverage*std_PMG0_m22,initSigmaCoverage*std_PMG0_m22,PMG_max_r)+meanX01(2);
    % PMG_deltaN = PMG_Nvec(2)-PMG_Nvec(1);
    % [PMG_Emat, PMG_Nmat] = meshgrid(PMG_Evec, PMG_Nvec);
    % ---
    % grid_delta = PMG_deltaE*PMG_deltaN;
    % grid_rbpmf = [PMG_Emat(:), PMG_Nmat(:)]';

    [grid_rbpmf, grid_delta, gridDimOldRb, xOldRb, PpoldRb, ~] = gridCreation(meanX0(1:2),varX0(1:2,1:2),sFactor,2,Npa(1:2));
    PMG_Emat = grid_rbpmf(1,:);
    PMG_Nmat = grid_rbpmf(2,:);
    grid_delta = prod(grid_delta);
    % - nonlinear part - PM PDF creation
    PMGpred = 1/((2*pi)*sqrt(Cov_PMG0_m(1,1))*sqrt(Cov_PMG0_m(2,2))) ...
        *exp(-0.5 *(((PMG_Emat(:)-meanX0(1)).^2)/Cov_PMG0_m(1,1)+((PMG_Nmat(:)-meanX0(2)).^2)/Cov_PMG0_m(2,2) ))';
    % - linear part
    dimLinState = nl;
    xpKF  = zeros(dimLinState, Nrb);
    for i=1:Nrb
        xpKF(:,i) = meanX02;
    end
    PpKF   = zeros(dimLinState, dimLinState, Nrb);
    for i=1:Nrb
        PpKF(:,:,i) = Pinit22;
    end
    % --
    % - global moments computation
    % -- linear part
    [xp_LIN, Pp_LIN] = GMmoments(PMGpred*grid_delta,xpKF,PpKF);
    % -- nonlinear part
    xp_NLIN = pdfMean(PMGpred,grid_rbpmf,grid_delta,nn,Nrb);
    Pp_NLIN = pdfCov(PMGpred,grid_rbpmf,xp_NLIN,grid_delta,nn,Nrb);
    % -- cross-covariance matrix
    crossLNp = crossCov(PMGpred, grid_rbpmf, xp_NLIN, xpKF, xp_LIN, grid_delta, nl, nn, Nrb);
    % --- moments
    xp_rbpmf(:,1) = [xp_NLIN; xp_LIN];
    Pp_rbpmf(:,:,1) = [Pp_NLIN, crossLNp'; crossLNp, Pp_LIN];
    % --
    %upmf = zeros(nn,1);
    % --
    % -- lin. state initialisation
    %zstar = zeros(nn,Nrb,Nrb);
    xpKFaux = zeros(nl,Nrb);
    PpKFaux = zeros(nl,nl,Nrb);
    % -- weights computation and Gaussian pdfs of linear state computation for
    %   the new grid
    % wp = zeros(Nrb,Nrb);       % 'transition' matrix of weights between old (filtering) and new (predictive) grids
    wpSummed = zeros(1,Nrb); % weights for new grid
    % xpKFnewGrid = zeros(nl,Nrb);
    % PpKFnewGrid = zeros(nl,nl,Nrb);
    % -- 'full' convolution (might be optimised)
    % --- relation between old and new grid (grid for state at time k and grid for time k+1)
    % --
    Qaux = Q11;% + 0.5*gridNew_delta*eye(nn); % overbounding
    invQaux = inv(Qaux);
    % --


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
            [predMean2(:,k+1), predVar2(:,:,k+1)] = momentGbF(predGrid, predDensityProb, predGridDelta(:,k+1));

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
            [predMean3(:,k+1), predVar3(:,:,k+1)] = momentGbF(predGrid2, predDensityProb2, predGridDelta2(:,k+1));

            tocPMF2(k) = toc;
        end

        %% RBPMF
        if 1%runFlags.LGbF_Spectral
            tic

            % measurement update
            % - preparation of measurement function
            % grid_rbpmf = [PMG_Emat(:), PMG_Nmat(:)]';
            hpmf = h1(grid_rbpmf);
            Hkf = [H2];
            % - PMF part
            for i=1:Nrb
                % --- innovation (respected linear term)
                innovation(:,i) = z(:,k) - (hpmf(:,i) + Hkf*xpKF(:,i));
                % --- likelihood
                likelihood(:,i) = mvnpdf1(innovation(:,i), zeros(nz,1), R, invR, nz);
            end
            % --- filtering pdf unnormalized
            PMGfilt = likelihood.*PMGpred;
            % --- normalisation
            alpha = sum(PMGfilt)*grid_delta;
            % --- filtering pdf
            PMGfilt = PMGfilt/alpha;
            % - KF part
            % -- filtering pdfs computation
            % --- KF for each point (a priori stat. of all points are the SAME)
            for i=1:Nrb
                %Kkf = KFsPred_var(:,:,i)*Hkf'*inv(Hkf*KFsPred_var(:,:,i)*Hkf'+R);
                Kkf = PpKF(:,:,i)*Hkf'/(Hkf*PpKF(:,:,i)*Hkf'+R);
                % xfKF(:,i) = xpKF(:,i) + Kkf*(z(:,k) - (hpmf(:,i) + Hkf*xpKF(:,i)));
                xfKF(:,i) = xpKF(:,i) + Kkf*innovation(:,i);
                PfKF(:,:,i) = PpKF(:,:,i) -Kkf*Hkf*PpKF(:,:,i);
            end
            % --- filtering moments calculation
            % --- filtering pdf moments
            xf_NLIN = pdfMean(PMGfilt,grid_rbpmf,grid_delta,nn,Nrb);
            Pf_NLIN = pdfCov(PMGfilt,grid_rbpmf,xf_NLIN,grid_delta,nn,Nrb);
            % --- KF filtering pdf global moments
            wfKF = PMGfilt*grid_delta;
            wfKF = wfKF/sum(wfKF);
            [xf_LIN, Pf_LIN] = GMmoments(wfKF,xfKF,PfKF);
            % --- filtering cross-covariance matrix (between linear and nonlinear 'states')
            crossLNf = crossCov(PMGfilt, grid_rbpmf, xf_NLIN, xfKF, xf_LIN, grid_delta, nl, nn, Nrb);
            % --- moments
            xf_rbpmf(:,k) = [xf_NLIN; xf_LIN];
            Pf_rbpmf(:,:,k) = [Pf_NLIN, crossLNf'; crossLNf, Pf_LIN];

      
            [PMGfilt,gridDimOldRb, grid_delta, grid_rbpmf, eigVectRb,xfKF,PfKF] = measPdfInterp2(PMGfilt, ...
                gridDimOldRb,xOldRb,PpoldRb,Q(1:2,1:2),sFactor,2,Npa(1:2),xf_rbpmf(1:2,k),Pf_rbpmf(1:2,1:2,k) ,F(1:2,1:2),xfKF,PfKF);
            gridNew_rbpmf = F(1:2,1:2)* grid_rbpmf + u(1:2,k);
            gridNew_delta = prod(F(1:2,1:2)*grid_delta);
            grid_delta = prod(grid_delta);
            xOldRb = F(1:2,1:2) * xf_rbpmf(1:2,k) + u(1:2,k);
            PpoldRb = F(1:2,1:2) * eigVectRb;

            % time update
            xpKFaux =  F22*xfKF;
            PpKFaux = F22^2*PfKF + Q22;

            halfGridInd = ceil(length(gridNew_rbpmf)/2); % Index of middle point of predictive grid
            distPoints = (gridNew_rbpmf(:,halfGridInd)'-(gridNew_rbpmf)'); % Distance of transformed grid points to the new grid points
            convKer = ((exp(sum(-0.5*distPoints*invQ(1:2,1:2).*distPoints,2)))/predDenDenomW2)';% Convolution kernel values
            convKer = reshape(convKer,Npa(1:2)); % Convolution kernel values in tensor format

            %--- precompute the "basis" = p_i * Δx
            basis = reshape(PMGfilt, Npa(1:2)) * prod(grid_delta);

            %--- do the three un-normalized convolutions
            W0 = convn( basis,                    convKer, 'same' );    % Σ p_i Δx ⋅ K_ij
            M0 = convn((reshape(xpKFaux,Npa(1:2)) .* basis), convKer,'same');  % Σ x_i p_i Δx ⋅ K_ij
            S0 = convn(((reshape(PpKFaux,Npa(1:2)) + reshape(xpKFaux,Npa(1:2)).^2) ...
                .* basis), convKer,'same');                  % Σ (Var_i+μ_i²) p_i Δx ⋅ K_ij

            %--- now normalize exactly as in your loop:
            normFactor = sum(W0(:)) * prod(gridNew_delta);     % = Σ_j Σ_i p_i Δx ⋅ K_ij  · Δy

            W  = W0 / normFactor;          % this is PMGpred(j)

            %— conditional mean and variance at each grid-point j:
            xpKF_j = M0 ./ W0;             % E[x | transition→y_j]
            varKF_j = S0 ./ W0 - xpKF_j.^2;

            %— reshape back into your desired output shapes
            xpKF  = reshape(xpKF_j,  1, []);         % same as in the loop
            PpKF  = reshape(varKF_j, 1, 1, []);      % same as the GMmoments outpu

            %--- Predictive pmf on new grid
            PMGpred = reshape(W,1,[]);

            % -- reallocation
            grid_rbpmf = gridNew_rbpmf;
            grid_delta = gridNew_delta;

            tocPMF3(k) = toc;
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
    if 1
        [rmsePMF3(:,mc), astdPMF3(:,mc), annesPMF3(mc)] = calcRes(x, xf_rbpmf, Pf_rbpmf, k);
        tocPMF3avg(mc) = mean(tocPMF3);
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
if 1
    annes_PMFout3Pom = sum(annesPMF3) / (nx * mc * (k + 1));
    rmsePMFout3      = mean(rmsePMF3, 2);
    astdPMFout3      = mean(astdPMF3, 2);
    tocPMFavg3out    = mean(tocPMF3avg);
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
if 1
    rowNames{end+1} = 'RBGbF';
    vals{end+1} = [rmsePMFout3(1), rmsePMFout3(2), rmsePMFout3(3), ...
        astdPMFout3(1), astdPMFout3(2), astdPMFout3(3), ...
        annes_PMFout3Pom, tocPMFavg3out];
end

T2 = cell2mat(vals');
T2 = array2table(T2, ...
    'VariableNames', {'RMSE_x1','RMSE_x2','RMSE_x3','ASTD_1','ASTD_2','ASTD_3','ANNES','TIME'}, ...
    'RowNames', rowNames);

disp(T2)

% Trajectory Plots
if runFlags.LGF_Standard, plot(x(1,:),x(2,:)); hold on; plot(measMean2(1,:),measMean2(2,:)); end
if runFlags.LGbF_Spectral, plot(measMean3(1,:),measMean3(2,:)); end
if 1, plot(xf_rbpmf(1,:),xf_rbpmf(2,:)); end
legend('true','LGbF','SLGbF','RBGbF')

figure
subplot(1,2,1)
plot(x(1,1:k-1)-measMean2(1,1:k-1))
hold on
plot(x(1,1:k-1)-measMean3(1,1:k-1))
plot(x(1,1:k-1)-xf_rbpmf(1,1:k-1))
legend('LGbF','SLGbF','RBGbF')
subplot(1,2,2)
plot(x(2,1:k-1)-measMean2(2,1:k-1))
hold on
plot(x(2,1:k-1)-measMean3(2,1:k-1))
plot(x(2,1:k-1)-xf_rbpmf(2,1:k-1))
legend('LGbF','SLGbF','RBGbF')
