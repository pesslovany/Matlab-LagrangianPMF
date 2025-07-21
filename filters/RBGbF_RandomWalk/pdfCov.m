function P = pdfCov(PMpdf,PMgrid,m,deltas,dim,N)

% P = zeros(dim,dim);
% for i=1:N
%     P = P + (PMgrid(:,i)-m)*(PMgrid(:,i)-m)'*PMpdf(i)*prod(deltas);
% end

aux = PMgrid-m;
aux_w = aux.*kron(ones(dim,1),PMpdf)*prod(deltas);

P = aux*aux_w';