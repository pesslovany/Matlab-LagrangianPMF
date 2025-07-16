function [mu, Sigma] = momentRBPF(measXn, measXl, measPl, measW, nxn, nxl, xnid, xlid)

    noPartRbpf = length(measW);
    
    mun = measXn*measW;
    mul = measXl*measW;
    mu = zeros(nxn + nxl, 1);
    mu(xnid) = mun;
    mu(xlid) = mul;
    
    Sigman = zeros(nxn, nxn);
    Sigmal = zeros(nxl, nxl);
    Sigmanl = zeros(nxn, nxl);
    for ii = 1:noPartRbpf
        Sigman = Sigman + measW(ii)*(measXn(:, ii) - mun)*(measXn(:, ii) - mun)';
        Sigmal = Sigmal + measW(ii)*(measPl(:, :, ii) + (measXl(:, ii) - mul)*(measXl(:, ii) - mul)');
        Sigmanl = Sigmanl + measW(ii)*(measXn(:, ii) - mun)*(measXl(:, ii) - mul)';
    end
    Sigma = zeros(nxn + nxl, nxn + nxl);
    Sigma(xnid, xnid) = Sigman;
    Sigma(xlid, xlid) = Sigmal;
    Sigma(xnid, xlid) = Sigmanl;
    Sigma(xlid, xnid) = Sigmanl';

end