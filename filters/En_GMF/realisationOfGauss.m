function [ensemblePosterior] = realisationOfGauss(resampleInd,XkGSF,PkGSF)
%REALISATIONOFGAUSS Summary of this function goes here
%   Detailed explanation goes here
% One draw from each Gaussian N(XkGSF(:,inds), PkGSF(:,:,inds)) — fast path (SPD covariances)

ki  = numel(resampleInd);
nx = size(XkGSF,1);

L = zeros(nx,nx,ki,'like',XkGSF);
for i = 1:ki
    L(:,:,i) = chol(PkGSF(:,:,resampleInd(i)),'lower');     % assumes SPD (fastest)
end

ensemblePosterior = XkGSF(:,resampleInd) + reshape( ...
    pagemtimes(L, randn(nx,1,ki,'like',XkGSF)), nx, ki);
end