%%
%   Lagrangian based grid-based filter for nonlinear models, simulated dataset.
%   authors: Jakub Matousek, pesslovany@gmail.com
%            Jindrich Dunik, dunikj@kky.zcu.cz, University of West Bohemia
% See paper: TBD doi.

%% Parameters and system simulation
clc
clear variables
close all
format shortG

addpath(genpath(pwd)); % Add all files and subfolders in the current directory


MC = 1;

for mc = 1:1:MC

    modelChoose = 2; % Which model to initialize

    model = initModelNonLin(modelChoose); % Initialize model

    % Unpack model structure variables
    fields = fieldnames(model);
    for i = 1:numel(fields)
        eval([fields{i} ' = model.' fields{i} ';']);
    end

    clear model

    % PMF initialization and parameters

    % Initial grid
    [predGrid, predGridDelta, gridDimOld, gridCenter, gridRotation] = gridCreation(meanX0,varX0,sFactor,nx,Npa);

    % Initial PMD
    fixTerm = (predGrid-meanX0);
    denominator = sqrt((2*pi)^nx*det(varX0));
    predPdf = ((exp(sum(-0.5*fixTerm'/(varX0).*fixTerm',2)))/denominator); % Initial Gaussian Point mass density (PMD)

    % Auxiliary variables
    halfGrid = ceil(N/2); % Middle row of the TPM matrix index

    noPart = 20000; % no Particles
    ksiPrior = mvnrnd(meanX0,varX0,noPart)'; %particles prior


    % For cycle over time stes
    for k = 1:1:kf
        disp(['Step:', num2str(k), '/', num2str(kf-1)])

        % Grid-based Filter
        tic

        % Measurement update
        [filtPdf] = gbfMeasUpdate(predGrid,nz,k,z(:,k),V,predPdf,predGridDelta,hfunct); % Measurement update

        % Filtering mean and var
        filtMeanPMF(:,k) = predGrid*filtPdf*prod(predGridDelta); %#ok<*SAGROW> % Measurement update mean
        chip_ = (predGrid-filtMeanPMF(:,k));
        chip_w = chip_.*repmat(filtPdf',nx,1);
        filtVarPMF(:,:,k) = chip_w*chip_' * prod(predGridDelta); % Measurement update variance
        filtGrid = predGrid;

        % Meas PDF interp
        [filtPdf, gridDimOld, GridDelta, measGridNew, gridRotation, gridCenter] = measPdfInterpNonLin(filtPdf,...
            gridDimOld,gridCenter,gridRotation,Q,sFactor,nx,Npa,filtMeanPMF(:,k),filtVarPMF(:,:,k),f,f_inv,jacobian_f,u(:,k));

        % Time Update
        [predPdf,predGrid,predGridDelta] = gbfTimeUpdateFFTnonLin(f,filtPdf,measGridNew,GridDelta,Npa,nx,u(:,k),Q);
        tocPMF(k) = toc; % Time evaluation

        %%
        tic
        predThrMeasEq = hfunct(ksiPrior,zeros(nz,1),k); %Prediction density grid through measurement EQ
        pom = z(:,k)'-predThrMeasEq'; %Measurement - measurementEQ(Grid)
        w = pdf(V.pdf,pom);
        w = w/sum(w);
        xEst(:,k) = ksiPrior*w; % mean
        varEst(:,k) = ksiPrior.^2*w - xEst(:,k).^2; %var

        % cumW = cumsum(w);
        % randomN = rand(1,noPart);
        % I = binarySearch(cumW, randomN, 'first');
        % ksi = ksiPrior(:,I);

        % O(n) Resampling without binary search or find
        cumW = cumsum(w);
        thresholds = (rand + (0:noPart-1)) / noPart;
        resampleInd = zeros(1, noPart);
        cumInd = 1;
        for p = 1:noPart
            while thresholds(p) >= cumW(cumInd)
                cumInd = cumInd + 1;
            end
            resampleInd(p) = cumInd;
        end
        ksi = ksiPrior(:, resampleInd);

        Qval = Q(ksi(end,ceil(length(ksi)/2)));
        ksiPrior = f(ksi,u(:,k)) + mvnrnd(zeros(1,nx), Qval, noPart)';

        tocPF(k) = toc;

    end

    rmsePMF(mc) = mean(sqrt(mean((x(:,1:kf-1)-filtMeanPMF(:,1:kf-1)).^2,2))); %#ok<*SAGROW>
    tocPMF(mc) = mean(tocPMF);

    rmsePF(mc) = mean(sqrt(mean((x(:,1:kf-1)-xEst(:,1:kf-1)).^2,2))); %#ok<*SAGROW>
    tocPF(mc) = mean(tocPF);

end

T = table( ...
    [mean(rmsePMF); mean(rmsePF)], ...
    [mean(tocPMF);  mean(tocPF)], ...
    'VariableNames', {'RMSE', 'Time'}, ...
    'RowNames', {'PMF', 'PF'});

disp(T)

figure
hold on
plot(filtMeanPMF(1,:),filtMeanPMF(2,:),"LineWidth",2,"Color","red")
plot(xEst(1,:),xEst(2,:),"LineWidth",2,"Color","blue")
plot(x(1,:),x(2,:),Color="black",Marker=".")
legend('PMF','PF','true state')