function [resampleInd] = systematicResampling(w,noPart)
%SYSTEMATICRESAMPLING Summary of this function goes here
%   Detailed explanation goes here
% O(n) Resampling without binary search or find
cumW = cumsum(w);
thresholds = (rand + (0:noPart-1)) / noPart;
resampleInd = zeros(1, noPart);
cumInd = 1;
for p = 1:noPart
    while thresholds(p) >= cumW(cumInd)
        cumInd = cumInd + 1;
    end
    resampleInd(p) = cumInd;
end

end