function [meanVec, diagCov] = momentsCPD(tensorCPD, grid, gridDelta)
% Compute mean and diagonal covariance of a CPD over a structured grid
% INPUTS:
%   tensorCPD  - CPD tensor (ktensor)
%   grid       - 1xD cell array with per-dimension coordinate vectors
%   gridDelta  - Dx1 vector of grid step sizes
% OUTPUTS:
%   meanVec    - Dx1 mean vector
%   diagCov    - Dx1 vector (diagonal of covariance)

R = length(tensorCPD.lambda);
D = ndims(tensorCPD);
meanVec = zeros(D,1);
diagCov = zeros(D,1);

% Precompute sum over factors
factorSums = cellfun(@(U) sum(U,1), tensorCPD.u, 'UniformOutput', false);

volume = prod(gridDelta);  % grid neighborhood size

for d = 1:D
    xd = grid{d}(:);         % grid coordinates
    Ud = tensorCPD.u{d};     % loading vectors for dim d

    % Mean value
    Ex = xd' * Ud;           % 1 x R
    prodRest = ones(1, R);
    for n = [1:d-1, d+1:D]
        prodRest = prodRest .* factorSums{n};
    end
    meanVec(d) = volume * sum(tensorCPD.lambda' .* Ex .* prodRest);

    % E[x_d^2]
    Ex2 = (xd.^2)' * Ud;     % 1 x R
    Ex2_val = volume * sum(tensorCPD.lambda' .* Ex2 .* prodRest);

    % Diagonal of covariance matrix
    diagCov(d) = Ex2_val - meanVec(d)^2;
end

end
