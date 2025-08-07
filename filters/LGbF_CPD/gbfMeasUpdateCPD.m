function [measPdf] = gbfMeasUpdateCPD(predGrid, z, noiseV, predDensityProb, predGridStep, hfunct, mapping, Npa, Rdef)
%pmfMeasMix - measurement update matlab pdf class noise
%INPUTS:
% predGrid - predictive grid
% nz - dimension of measurement
% k - time step
% z - measurement
% noiseV - noiseV.pdf is a measurement noise by dimensions
% predDensityProb - predictive pdf in CPD
% predGridStep - predictive grid step
% hfunct - measurement equation by dimension
% mapping - which state dimmensions map to which measurements
% Npa - number of points per axis
%OUTPUTS:
% measPdf - measurement pdf

measPdf = predDensityProb;  % start with prior CPD

% Over each measurement
for indMeas = 1:length(mapping.indexed)

    % Extract relevant grid axes and form grid for dependent state
    % variables
    gridDims = mapping.Combs{mapping.indexed(indMeas)};
    gridMeas = combvec(predGrid(gridDims));

    % Calculate likelihood one measurement
    predGridTrsf = hfunct{indMeas}(gridMeas);
    predGridTrsf(1, isnan(predGridTrsf(1,:))) = inf;  % clip extrapolations
    inov = z(indMeas)' - predGridTrsf';
    likelihood = reshape(pdf(noiseV(indMeas).pdf, inov), Npa(gridDims));

    % Decompose likelihood
    switch length(gridDims)
        % If vector it is directly CPD
        case 1
            % likelihood must be column vector
            likeCPD = ktensor(1, likelihood);
        % If matrix SVD decompose
        case 2
            % Decompose likelihood into CPD format
            [U, S, V] = svd(likelihood, "econ", "vector");
            cumEnergy = cumsum(S.^2) / sum(S.^2);
            r = binarySearch(cumEnergy, 0.9999);
            U = U(:,1:r); S = S(1:r); V = V(:,1:r);
            likeCPD = ktensor(S, {U, V});
        % If tensor CPD decomposition
        otherwise
            likeCPD = cp_als(likelihood, Rdef);
    end

    % Construct CPD cores with ones on non changing dimensions
    R = height(S);
    D = ndims(predDensityProb);
    factors = cell(1, D);
    for d = 1:D
        idx = find(gridDims == d);
        if ~isempty(idx)
            factors{d} = likeCPD.u{idx};
        else
            factors{d} = ones(size(predDensityProb.u{d},1), R);
        end
    end

    % Elementwise multiplication of likelihood with weights
    measPdf = hadamard_ktensor(measPdf, ktensor(S, factors));

    % Deflate
    if ncomponents(measPdf) > Rdef
        measPdf = cp_als(measPdf,Rdef);
    end

end

% Normalize
Zsum = sumTensor(measPdf);
measPdf.lambda = measPdf.lambda / (Zsum * prod(predGridStep));


end