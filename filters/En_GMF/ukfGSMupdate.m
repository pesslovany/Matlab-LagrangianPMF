function [XkGSF,PkGSF,wkGSF] = ukfGSMupdate(Xbark,wbark,nx,Ps,z,k,hfunct,R)
%UKFGSMUPDATE Summary of this function goes here
%   Detailed explanation goes here
[s,n] = size(Xbark);
nz = height(z);
alpha = 1e-3; beta = 2; kappa = 0;
lambda = alpha^2*(nx + kappa) - nx;
c = nx + lambda;
S  = chol(c*Ps,'lower');
Sig = [zeros(nx,1) S -S];
L  = size(Sig,2);
Xi = repmat(reshape(Xbark,[nx 1 n]),[1 L 1]) + repmat(reshape(Sig,[nx L 1]),[1 1 n]);
Yi_flat = hfunct(reshape(Xi,nx,[]), zeros(nz,1), k+1);
Yi_flat(isnan(Yi_flat)) = intmax('int32');
Yi = reshape(Yi_flat, nz, L, n);
wm = [lambda/c, repmat(1/(2*c),1,2*nx)];
wc = wm; wc(1) = wc(1) + (1 - alpha^2 + beta);
zhat = pagemtimes(Yi, reshape(wm.',[L 1 1]));
Yc = Yi - zhat; Xc = Xi - repmat(reshape(Xbark,[nx 1 n]),[1 L 1]);
Ycw = Yc .* reshape(wc,1,[],1);
Pzz = pagemtimes(Ycw, permute(Yc,[2 1 3])) + R;
Xcw = Xc .* reshape(wc,1,[],1);
Pxz = pagemtimes(Xcw, permute(Yc,[2 1 3]));
K     = pagemrdivide(Pxz, Pzz);
v     = z(:,k+1) - reshape(zhat,nz,n);
XkGSF = Xbark + reshape(pagemtimes(K, reshape(v,nz,1,n)), nx, n); % Mean vals
Kt    = permute(K,[2 1 3]);
PkGSF = Ps - pagemtimes(pagemtimes(K,Pzz), Kt); % Covars
vt    = permute(reshape(v,nz,1,n),[2 1 3]);
% detP  = reshape(prod(pageeig(Pzz,'vector'),1), n, 1);
a = reshape(Pzz(1,1,:),[],1);
b = reshape(Pzz(1,2,:),[],1);
c = reshape(Pzz(2,1,:),[],1);
d = reshape(Pzz(2,2,:),[],1);
logdetP = log(a.*d - b.*c);  % one imag nubmer
quad  = reshape(-0.5*pagemtimes(vt, pagemldivide(Pzz, reshape(v,nz,1,n))), n, 1);
wkGSF = log(wbark) - 0.5*logdetP + quad;
m     = max(wkGSF);
wkGSF = exp(wkGSF - (m + log(sum(exp(wkGSF - m)))));
wkGSF = wkGSF / sum(wkGSF); % Weights

end