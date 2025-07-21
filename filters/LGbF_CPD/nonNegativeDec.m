function [Ux,Sx,Vx] = nonNegativeDec(Zx,varargin)
%NN_DECOMP   Non-negative low-rank decomposition with automatic rank pick
%
%   [Ux,Sx,Vx,sx,cumEx,Kx,Ux_tr,Vx_tr] = nn_decomp(Zx)
%   automatically chooses Kx such that the first Kx singular values of Zx
%   (from an initial SVD) capture ≥99.999% of the energy, then runs
%       [W,H] = nnmf(Zx,Kx)
%   and returns
%     Ux      = W (MxKx nonnegative basis)
%     Vx      = H' (NxKx nonnegative coefficients)
%     Sx      = eye(Kx)  (so that Ux*Sx*Vx' = W*(H))
%     sx      = diag(S_svd)  (for diagnostics only)
%     cumEx   = cumulative energy fraction of the svd
%     Kx      = chosen rank
%     Ux_tr   = Ux
%     Vx_tr   = Vx
%
%   [...] = nn_decomp(Zx,ENERGY_THRESH) lets you override the 0.99999
%   threshold (default = 0.99999).
%
%   Example:
%     [Ux,Sx,Vx]    = deal([]);  % just to show names
%     [Ux,Sx,Vx,sx,cumEx,Kx,Ux_tr,Vx_tr] = nn_decomp(Zx,0.9999);
%

% parse optional threshold
if isempty(varargin)
    energyThresh = 0.9999;
else
    energyThresh = varargin{1};
end

% 1) quick SVD to pick rank
[~,S_svd,~] = svd(Zx,'econ');
sx   = diag(S_svd);
cumEx = cumsum(sx.^2)/sum(sx.^2);
Kx   = find(cumEx>=energyThresh,1,'first');

% 2) run non-negative MF of that rank
opts = statset('maxiter',1000,'display','off');
[W,H] = nnmf(Zx, Kx, 'options', opts);

% 3) package outputs with same names
Ux     = W;           % M×K
Vx     = H';          % N×K
Sx     = eye(Kx);     % K×K


end