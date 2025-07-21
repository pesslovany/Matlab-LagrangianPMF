function [cpd] = gbfDiffusionCPDSpect(predDensityCPD, Q, predGridDeltaCPD, ...
    axisp, Npa, nx, k, Rdef)
%GBFDIFFUSIONCPD diffusion using spectral methods in CPD format

% 1) Precompute lengths L and diagonal-only Uhat
L = cellfun(@(v) v(end)-v(1), axisp);
L = L + predGridDeltaCPD(:,k+1);
Uhat = cell(nx,1);
for d = 1:nx
    % build kvec and attenuation kernel
    kk = fftshift(floor(-(Npa(d))/2:(Npa(d))/2-1)') * (2*pi / L(d));
    Kd = exp( - ((Q(d,d)/2) .* (kk.^2)));

    % get original factor and apply Dirichlet BC
    U = predDensityCPD.U{d};
    U([1 end],:) = 0;

    % FFT all columns at once, apply Kd, then IFFT
    Uhat{d} = real( ifft( bsxfun(@times, Kd, fft(U, Npa(d),1)), ...
        Npa(d), 1 ) );
end

% assemble diag-only result
cpd = ktensor(predDensityCPD.lambda, Uhat);
cpd_diag = cpd;

% 2) Precompute FFTs of each diag factor for cross-terms
FFTdiag = cell(nx,1);
for d = 1:nx
    FFTdiag{d} = fft(cpd.U{d}, Npa(d), 1);  % [NpaCPD(d) x R]
end

% 3) Loop over non-zero off-diagonals only
pairs = find(triu(Q,1));    % linear indices of i<j where Q(i,j)~=0
[Is, Js] = ind2sub([nx nx], pairs);

for p = 1:numel(Is)
    i = Is(p);  j = Js(p);
    qij = Q(i,j);
    if qij == 0, continue; end

    % build 2D spectral grid and do economical SVD
    ki = fftshift(floor(-(Npa(i))/2:(Npa(i))/2-1)') * (2*pi / L(i));
    kj = fftshift(floor(-(Npa(j))/2:(Npa(j))/2-1)') * (2*pi / L(j));
    [Ki,Kj] = ndgrid(ki,kj);
    Aij = exp(-(qij .* (Ki .* Kj)));

    % [Ui,Si,Vj] = svd(Aij,'econ');
    % si     = diag(Si);
    % cumEx  = cumsum(si.^2)/sum(si.^2);
    % rsvd     = find(cumEx>=0.99999,1,'first');
    % Ui = Ui(:,1:rsvd);
    % Vj = Vj(:,1:rsvd);
    % sigma = diag(Si(1:rsvd,1:rsvd));

    [Ui,sigma,Vj] = nonNegativeDec(Aij);
    sigma = diag(sigma);

    % Now vectorize the inner r-loop:
    % Instead multiply per ℓ:
    for ell = 1:length(sigma)
        h_i = Ui(:,ell);
        h_j = Vj(:,ell);

        % apply filter to all R components at once:
        Ui_filt = real( ifft( bsxfun(@times, h_i,  FFTdiag{i} ), Npa(i),1 ) );
        Uj_filt = real( ifft( bsxfun(@times, h_j,  FFTdiag{j} ), Npa(j),1 ) );
        Ui_filt([1 end],:) = 0;
        Uj_filt([1 end],:) = 0;

        % build a temporary cell array with *all* updated modes in one go
        % by copying cpd_diag.U and replacing mode i and j
        tempFactors = cpd_diag.U;
        tempFactors{i} = Ui_filt;
        tempFactors{j} = Uj_filt;

        % new weights: vector of length R
        newLambdas = sigma(ell) * predDensityCPD.lambda;

        % add all rank-1 terms at once:
        cpd = cpd + ktensor(newLambdas, tempFactors);
    end
end

if ncomponents(cpd) > Rdef
    cpd = cp_als(cpd, Rdef);
end

Zsum       = sumTensor(cpd);
cpd.lambda = cpd.lambda / (Zsum * prod(predGridDeltaCPD(:,k+1)));

end