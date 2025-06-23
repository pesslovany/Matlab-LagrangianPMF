%%
%   Lagrangian grid based filter based on FFT and spectral differentiation
%   for real TAN dataset.
%   author: pesslovany@gmail.com
%   See paper (in review will be changed)
%% Parameters and system simulation
clc
clear variables
close all
format shortG

addpath(genpath(pwd)); % Add all files and subfolders in the current (one level up) directory

modelChoose = 4; % choose model 3D - 3, 4D with - 4

load('data.mat') % map of terrain
% Map interpolant for measurement
vysky = griddedInterpolant(souradniceX',souradniceY',souradniceZ',"linear","none");
% Time steps
timeSteps = 1:1:length(souradniceGNSS);
% Last time
endTime = length(timeSteps);


% Number of Monte Carlo simulations
MC = 50;
for mc = 1:1:MC

    model = initModel(modelChoose, souradniceGNSS, hB, vysky);

    % Unpack model structure variables
    fields = fieldnames(model);
    for i = 1:numel(fields)
        eval([fields{i} ' = model.' fields{i} ';']);
    end

    clear model

    dtSpec = 0.01;
    sFactor = 6;

    % PMF init and param
    % Initial grid
    [predGrid, predGridDelta, gridDimOld, xOld, Ppold, ~] = gridCreation(meanX0,varX0,sFactor,nx,Npa);
    [predGrid2, predGridDelta2, gridDimOld2, xOld2, Ppold2, ~] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial PMD
    pom = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predDensityProb = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)
    predDensityProb2 = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)
    predPdf3 = ((exp(sum(-0.5*pom'/(varX0).*pom',2)))/denominator); % Initial Gaussian Point mass density (PMD)
    % Auxiliary variables
    predDenDenomW = sqrt((2*pi)^nx*det(Q)); % Denominator for convolution in predictive step
    halfGrid = ceil(N/2); % Middle row of the TPM matrix index
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
        measGrid = predGrid;
        [measPdf,gridDimOld, GridDelta(:,k), measGridNew, eigVect] = measPdfInterp(measPdf,...
            gridDimOld,xOld,Ppold,Q,sFactor,nx,Npa,measMean2(:,k),measVar2(:,:,k),F);
        xOld = F*measMean2(:,k) + u(:,k);
        Ppold = F*(eigVect);
        % Time Update
        [predDensityProb,predGrid,predGridDelta] = gbfTimeUpdateFFT(F,measPdf,measGridNew,GridDelta,k,Npa,invQ,predDenDenomW,nx,u(:,k));
        tocPMF(k) = toc; % Time evaluation

        %% Grid-based Filter - cont
        tic
        % Measurement update
        [measPdf2] = gbfMeasUpdate(predGrid2,nz,k,z(:,k),V,predDensityProb2,predGridDelta2(:,k),hfunct); % Measurement update
        % Filtering mean and var
        measMean3(:,k) = predGrid2*measPdf2*prod(predGridDelta2(:,k)); % Measurement update mean
        chip_ = (predGrid2-measMean3(:,k));
        chip_w = chip_.*repmat(measPdf2',nx,1);
        measVar3(:,:,k) = chip_w*chip_' * prod(predGridDelta2(:,k)); % Measurement update variance
        measGrid2 = predGrid2;

        gridCenter = xOld2;
        gridRotation = Ppold2;
        gridDim = gridDimOld2;
        predVarPom = F*measVar3(:,:,k)*F' + Q;
        xOld2 = F*measMean3(:,k) + u(:,k);
        [predGridAdvect,  GridDelta2(:,k+1), gridDimOld2, ~, Ppold2, gridBound] = gridCreation(xOld2,predVarPom,sFactor,nx,Npa); % create the grid to interp on
        measGridNew2 = inv(F)*(predGridAdvect - u(:,k));
        GridDelta2(:,k) = inv(F)* GridDelta2(:,k+1);

        % Interpolation
        Fint = griddedInterpolant(gridDim,reshape(measPdf2,Npa),"linear","none"); % interpolant
        inerpOn = inv(gridRotation)*(measGridNew2 - gridCenter); % Grid to inter on transformed to rotated space
        measPdf2 = Fint(inerpOn'); % Interpolation
        measPdf2(isnan(measPdf2)) = 0; % Zeros for extrapolation, otherwise artifacts would appear

        rotQ = Ppold2'*Q*Ppold2;
        % Time Update
        [predDensityProb2,predGrid2,predGridDelta2] = gbfTimeUpdateSpect(F,measPdf2,measGridNew2,GridDelta2,k,Npa,rotQ,u(:,k),dtSpec,gridBound);
        % If something unexpected happend switch to discrete diffusion for
        % one step as it is more stable
        if abs(min(predDensityProb2)) > max(predDensityProb2)/100
            [predDensityProb2,predGrid2,predGridDelta2] = gbfTimeUpdateFFT(F,measPdf2,measGridNew2,GridDelta2,k,Npa,invQ,predDenDenomW,nx,u(:,k));
        end
        predDensityProb(predDensityProb<0)=0;
        predDensityProb2 = predDensityProb2./(sum(predDensityProb2)*prod(GridDelta2(:,k+1)))'; % Normalizaton (theoretically not needed)

        tocPMF2(k) = toc; % Time evaluation

    end
    % Evaluation
    rmsePMF(:,mc) = sqrt(mean((x(:,1:k-1)-measMean2(:,1:k-1)).^2,2)); %#ok<*SAGROW>
    astdPMF11(mc) = sqrt(mean(measVar2(1,1,1:k-1)));
    astdPMF22(mc) = sqrt(mean(measVar2(2,2,1:k-1)));
    astdPMF33(mc) = sqrt(mean(measVar2(3,3,1:k-1)));
    annes_PMF = 0;
    for indAn = 1:1:k-1
        annes_PMF = annes_PMF + ((x(:,indAn)-measMean2(:,indAn)).*(1./(diag(measVar2(:,:,indAn)))))'*(x(:,indAn)-measMean2(:,indAn));
    end
    annes_PMFout(mc) = annes_PMF;

    % Evaluation
    rmsePMF2(:,mc) = sqrt(mean((x(:,1:k-1)-measMean3(:,1:k-1)).^2,2)); %#ok<*SAGROW>
    astdPMF211(mc) = sqrt(mean(measVar3(1,1,1:k-1)));
    astdPMF222(mc) = sqrt(mean(measVar3(2,2,1:k-1)));
    astdPMF233(mc) = sqrt(mean(measVar3(3,3,1:k-1)));
    annes_PMF2 = 0;
    for indAn = 1:1:k-1
        annes_PMF2 = annes_PMF2 + ((x(:,indAn)-measMean3(:,indAn)).*(1./(diag(measVar3(:,:,indAn)))))'*(x(:,indAn)-measMean3(:,indAn));
    end
    annes_PMFout2(mc) = annes_PMF2;

    tocPMF2avg(mc) = mean(tocPMF2);

    tocPMFavg(mc) = mean(tocPMF);
end


annes_PMFout = sum(annes_PMFout)/(nx*mc*(k+1));
annes_PMFout2 = sum(annes_PMFout2)/(nx*mc*(k+1));

rmsePMFout = mean(rmsePMF,2);
rmsePMFout2 = mean(rmsePMF2,2);
tocPMFavg2out = mean(tocPMF2avg,2);
tocPMFavgOut = mean(tocPMFavg,2);

T2 = table([ rmsePMFout(1) rmsePMFout2(1)]',...
    [ rmsePMFout(2) rmsePMFout2(2)]',...
    [ rmsePMFout(3) rmsePMFout2(3)]',...
    [ mean(astdPMF11) mean(astdPMF211)]',...
    [ mean(astdPMF22) mean(astdPMF222)]',...
    [ mean(astdPMF33) mean(astdPMF233)]',...
    [ mean(annes_PMFout) mean(annes_PMFout2)]',...
    [ mean(tocPMFavgOut) mean(tocPMFavg2out)]',...
    'VariableNames',{'RMSE x1','RMSE x2','RMSE x3','ASTD 1','ASTD 2','ASTD 3','ANNES','TIME'},'RowName',...
    {'LGbF','Spect LGbF'}) %#ok<NOPTS>

plot(x(1,:),x(2,:))
hold on
plot(measMean2(1,:),measMean2(2,:))
plot(measMean3(1,:),measMean3(2,:))
