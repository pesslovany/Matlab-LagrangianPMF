function [predPdf,predGrid,gridStepOut,convKerTens] = gbfTimeUpdateFFTnonLin(f,measPdf,measGrid,gridStepIn,Npa,nx,u,Q)
%PMFUPDATEFFT time update in the Lagrangian formulation
%INPUTS:
% f - model dynamics matrix
% measPdf - measurement pdf
% measGrid - measurement grid
% gridStep - measurement grid step per dimension (and all old grid steps)
% k - time step
% Npa - number of points per axis
% nx - dimension of the state
% u - input to the model dynamics
%OUTPUTS:
% predPdf - predictive pdf
% predGrid - predictive grid
% gridStep - all grids steps with added predictive grid step
% convKerTens - convolution kernel in frequency space

% Pred Grid
predGrid = f(measGrid,u); % Predictive grid
gridStepOut = abs(predGrid(:,2) - predGrid(:,1)); % Predictive grid step size

measPdfDotDeltas = (measPdf.*gridStepIn'); % measurement PDF * measurement PDF step size
mesPdfDotDeltasTens = reshape(measPdfDotDeltas,Npa); % Into tensor space

halfGridInd = ceil(length(predGrid)/2); % Index of middle point of predictive grid
distPoints = (predGrid(:,halfGridInd)'-(predGrid)'); % Distance of transformed grid points to the new grid points

gridForQ = predGrid(end,halfGridInd);
Qval = Q(gridForQ); % Q is 5x5xN
% Compute determinant for each slice
predDenDenomW = sqrt((2*pi)^nx * det(Qval)); % 1xN
% Compute inverse of Q
invQ = inv(Qval); % 5x5xN
convKer = ((exp(sum(-0.5*distPoints*invQ.*distPoints,2)))/predDenDenomW)';% Convolution kernel values

convKerTens = reshape(convKer,Npa); % Convolution kernel values in tensor format

[predDensityProb2cub, convKerTens] = convFftN(mesPdfDotDeltasTens, convKerTens, Npa, 0); % FFT convolution

predPdf = reshape(predDensityProb2cub,length(predGrid),1); % back to vector for easier manipulation
predPdf = predPdf./(sum(predPdf)*prod(gridStepOut))'; % Normalizaton (theoretically not needed)
predPdf(predPdf<0) = 0; % FFT approach sometimes can produce very small negative values

end

