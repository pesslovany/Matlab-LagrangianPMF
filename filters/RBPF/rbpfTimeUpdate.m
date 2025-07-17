function [predXn, predXl, predPl, predW] = rbpfTimeUpdate(measXn, measXl, measPl, measW, F, Q, u, nxn, nxl, xnid, xlid)

noPartRbpf = length(measW);

Fnn = F(xnid, xnid);
Fnl = F(xnid, xlid);
Fln = F(xlid, xnid);
Fll = F(xlid, xlid);
un = u(xnid);
ul = u(xlid);
Qnn = Q(xnid, xnid);
Qnl = Q(xnid, xlid);
Qln = Q(xlid, xnid);
Qll = Q(xlid, xlid);

% PF
sqrtQ = chol(Q, "lower");
measXlpert = zeros(nxl, noPartRbpf);
for ii = 1:noPartRbpf
    measXlpert(:, ii) = chol(measPl(:, :, ii), "lower")*randn(nxl, 1);
end
noisen = sqrtQ(xnid, :)*randn(nxn + nxl, noPartRbpf);
FnlXlpert = Fnl*measXlpert;
predXn = Fnn*measXn + (Fnl*measXl + FnlXlpert) + un + noisen;
predW = measW;

% KF
predXl = zeros(size(measXl));
predPl = zeros(size(measPl));
Abar = Fll - Qln*inv(Qnn)*Fnl;
Qllbar = Qll - Qln*inv(Qnn)*Qnl;
pseInov = FnlXlpert + noisen;
pseZ = Fnl*measXl + pseInov;
for ii = 1:noPartRbpf
    N = Fnl*measPl(:, :, ii)*Fnl' + Qnn;
    L = Abar*measPl(:, :, ii)*Fnl'*N;
    predPl(:, :, ii) = Abar*measPl(:, :, ii)*Abar' + Qllbar - L*N*L';
    predXl(:, ii) = Abar*measXl(:, ii) + Qln*inv(Qnn)*pseZ(:, ii) + Fln*measXn(:, ii) + L*pseInov(:, ii) + ul;
end

end