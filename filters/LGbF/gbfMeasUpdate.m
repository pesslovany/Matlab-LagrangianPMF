function [measPdf] = gbfMeasUpdate(predGrid, nz, k, z, noiseV, predDensityProb, predGridStep, hfunct)
%pmfMeasMix - measurement update matlab pdf class noise
%INPUTS:
% predGrid - predictive grid
% nz - dimension of measurement
% k - time step
% z - measurement
% noiseV - noiseV.pdf is a measurement noise Matlab pdf
% predDensityProb - predictive pdf
% predGridStep - predictive grid step
% hfunct - measurement equation
%OUTPUTS:
% measPdf - measurement pdf

predGridTrsf = hfunct(predGrid,zeros(nz,1),k); %Prediction density grid through measurement EQ
predGridTrsf(1,isnan(predGridTrsf(1,:))) = inf; % To protect from map extrapolations
inov = z'-predGridTrsf'; %Measurement - measurementEQ(Grid)
measPdfNoNorm =  pdf(noiseV.pdf,inov).*predDensityProb; % Filtration density unnormalized
measPdf = measPdfNoNorm/(prod(predGridStep)*sum(measPdfNoNorm)); %Normalization

end
