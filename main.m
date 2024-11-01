%%
%   Lagrangian based grid-based filter for real TAN data.
%   authors: Jakub Matousek, pesslovany@gmail.com
%            Jindrich Dunik, dunikj@kky.zcu.cz, University of West Bohemia
% See paper: 10.1109/MSP.2024.3489969.

%% Parameters and system simulation
clc
clear variables
close all
format shortG

addpath(genpath(pwd)); % Add all files and subfolders in the current (one level up) directory

modelChoose = 3; % choose model 3D - 3, 4D with 2D measurement - 4

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

        % Grid-based Filter
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

    tocPMFavg(:,mc) = mean(tocPMF);

end

% Normalize outputs
annes_PMFout = annes_PMFout / (nx * k);
rmsePMFout = mean(rmsePMF, 2);
tocPMFavgOut = mean(tocPMFavg, 2);

