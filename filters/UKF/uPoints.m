function [chi, w] = uPoints(x, P)
    nx = size(x,1);
    S = chol(P, 'lower');
    decomp = sqrt(nx) * S;
    chi = repmat(x, 1, 2*nx+1);
    chi(:, 2:nx+1)       = x + decomp;
    chi(:, nx+2:end)     = x - decomp;
    w = ones(1, 2*nx+1) / (2*nx + 1);
end