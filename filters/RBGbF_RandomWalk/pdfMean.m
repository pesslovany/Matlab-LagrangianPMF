function m = pdfMean(PMpdf,PMgrid,deltas,dim,N)

% m = zeros(dim,1);
% for i=1:N
%     m = m + PMgrid(:,i)*PMpdf(i)*prod(deltas);
% end

m = PMgrid*PMpdf'*prod(deltas);