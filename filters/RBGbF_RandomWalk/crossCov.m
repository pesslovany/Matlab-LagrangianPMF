function P = crossCov(PMpdf,PMgrid,PMmean,KFestimates,KFmean,deltas,dimL,dimN,N)

% P = zeros(dimL,dimN);
% for i=1:N
%     P = P + (KFestimates(:,i))*(PMgrid(:,i))'*PMpdf(i)*prod(deltas);
% end
% P = P - KFmean*PMmean';

aux_w = PMgrid.*kron(ones(dimL,1),PMpdf)*prod(deltas);
P = KFestimates*aux_w' - KFmean*PMmean';
