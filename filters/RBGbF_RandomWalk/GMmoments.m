function [GMmean,GMvar] = GMmoments(weights, Means, Vars)

GMmean = weights*Means';
GMvar = weights*(squeeze(Vars)+(Means'-GMmean).^2);