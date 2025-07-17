function [predDensityProb,predGrid,GridDelta] = gbfTimeUpdateSpect(F,measPdf,measGridNew,GridDelta,k,Npa,Q,u,gridDimOld)
% LGbF time update for the spectral method

% Pred Grid
predGrid = F*measGridNew + u; % Predictive grid

% ULTRA FAST PMF
filtDenDOTprodDeltas = (measPdf*prod(GridDelta(:,k))); % measurement PDF * measurement PDF step size
filtDenDOTprodDeltasCub = reshape(filtDenDOTprodDeltas,Npa); % Into physical space

filtDenDOTprodDeltasCub = setEdgesToZeros(filtDenDOTprodDeltasCub);

dims = numel(gridDimOld);

L = 2*gridDimOld + GridDelta(:,k+1); % Grid size

kInd = cell(1, dims);
gridSize = zeros(1, dims);

for d = 1:dims
    gridSize(d) = Npa(d);
    kInd{d} = 2 * pi / L(d) * fftshift(floor(-(Npa(d)) / 2 : (Npa(d)) / 2 - 1)).';
end

% Initialize coefficient with ones of the final tensor shape
coeficienty = zeros(gridSize);

% Add diagonal terms (squared terms for each dimension)
for d = 1:dims
    shape = ones(1, dims);
    shape(d) = gridSize(d);
    coeficienty = coeficienty + reshape(kInd{d}.^2 * (Q(d, d) / 2), shape);
end

% Add cross-terms (products of different dimensions)
for d1 = 1:dims
    for d2 = d1+1:dims
            shape1 = ones(1, dims);
            shape2 = ones(1, dims);
            shape1(d1) = gridSize(d1);
            shape2(d2) = gridSize(d2);
            coeficienty = coeficienty + Q(d1,d2)*reshape(kInd{d1}, shape1) ...
                .* reshape(kInd{d2}, shape2);
    end
end

coeficienty = exp(-coeficienty);

% dims = 1:1:nx;
u = filtDenDOTprodDeltasCub;
u = fftn(u);
u = (coeficienty).*u;
predDensityProb2cub = real(ifftn(u)); % realna cast

predDensityProb = reshape(predDensityProb2cub,length(predGrid),1); % back to computational space

end

