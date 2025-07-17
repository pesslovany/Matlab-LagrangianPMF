function [measXn, measXl, measPl, measW] = rbpfMeasUpdate(predXn, predXl, predPl, predW, nx, nxl, nz, xnid, xlid, k, z, R, hfunct, hJacobLfunct, essThrdRbpf)
% Assumptions
% - Measurement noise ~ unimodal Gaussian

noPartRbpf = length(predW);
% PF
S = zeros(nz, nz, noPartRbpf);
invS = zeros(nz, nz, noPartRbpf);
predX = zeros(nx, noPartRbpf);
predX(xnid, :) = predXn;
predX(xlid, :) = predXl;
predhn = hfunct(predX, zeros(nz,1), k);
jacobL = zeros(nz, nxl, noPartRbpf);
inov = z - predhn;
measW = zeros(size(predW));
for ii = 1:noPartRbpf
    jacobL(:, :, ii) = hJacobLfunct(predX(:, ii));
    HPH = jacobL(:, :, ii)*predPl(:, :, ii)*jacobL(:, :, ii)';
    S(:,:,ii) = HPH + R;
    invS(:,:,ii) = inv(S(:,:,ii));
    
    measW(ii) = exp(-(1/2)*inov(:,ii)'*invS(:,:,ii)*inov(:,ii))*predW(ii);
end
measW = measW/sum(measW);
measXn = predXn;
% Resampling
if(1/sum(measW.^2) < essThrdRbpf)
    ind = sysresample(measW);
    measW = 1/noPartRbpf*ones(noPartRbpf,1);
    measXn = measXn(:,ind);
    predXl = predXl(:,ind);
    predPl = predPl(:,:,ind);
    S = S(:,:,ind);
    invS = invS(:,:,ind);
    inov = inov(:,ind);
end

% KF
measXl = zeros(size(predXl));
measPl = zeros(size(predPl));
for ii = 1:noPartRbpf
    K = (predPl(:,:,ii)*jacobL(:, :, ii)')*invS(:,:,ii);
    measPl(:,:,ii) = predPl(:,:,ii) - K*S(:,:,ii)*K';
    measPl(:,:,ii) = (measPl(:,:,ii) + measPl(:,:,ii)')*0.5;
    measXl(:,ii) = predXl(:,ii) + K*inov(:,ii);
end

end