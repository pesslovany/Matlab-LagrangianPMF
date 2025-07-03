function [mu, Sigma] = momentGbF(gridPoints, pdfVals, gridDelta)
%MOMENTGBF  Compute mean and covariance of a grid-based filter (GbF) PDF
%
%   [mu, Sigma] = momentGbF(gridPoints, pdfVals, gridDelta)
%
%   INPUTS
%     gridPoints : [nx × N] matrix – each column is a state point
%     pdfVals    : [1 × N]  vector – weights at the grid points
%     gridDelta  : [nx × 1] (or scalar) – step size per dimension
%
%   OUTPUTS
%     mu    : [nx × 1] mean vector
%     Sigma : [nx × nx] covariance matrix
%
%   Assumes ∑(pdfVals) ≈ 1 / prod(gridDelta).
%
% -------------------------------------------------------------------------

    dv = prod(gridDelta);          % volume element
    mu = gridPoints * pdfVals(:) * dv;

    delta = gridPoints - mu;       % centred coordinates
    Sigma = (delta .* pdfVals(:)') * delta.' * dv;
end
