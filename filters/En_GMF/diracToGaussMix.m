function [Ps] = diracToGaussMix(Xbark,wbark,Q)
%DIRACTOGAUSSMIX Summary of this function goes here
%   Detailed explanation goes here

[s,n] = size(Xbark);
xbark = Xbark*(wbark);
chip_ = (Xbark-xbark);
% Compute Ps using optimal weighting
alpha = 1; % adjust alpha if needed
Ps = alpha * (4 / ((n) * (s + 2)))^(2 / (s + 4)) * (chip_ .* repmat(wbark', s, 1)) * chip_' + Q;
Ps = (Ps + Ps.') / 2;

end