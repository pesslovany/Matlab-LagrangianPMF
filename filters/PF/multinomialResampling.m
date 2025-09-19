function [I] = multinomialResampling(w,noPart)
%MULTINOMIALRESAMPLING Summary of this function goes here
%   Detailed explanation goes here
% Resampling
cumW = cumsum(w);
randomN = rand(1,noPart);
I = binarySearch(cumW, randomN, 'first');
end