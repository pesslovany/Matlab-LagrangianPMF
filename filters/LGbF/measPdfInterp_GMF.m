function [measPdfNew,gridDimNew, gridStepNew, measGridNew, eigVectNew] = measPdfInterp_GMF(measPdf,gridDim,center,eigVect,Q,sFactor,nx,Npa,measMean,measVar,F,meanNew,varNew)
% Interpolates the measurement pdf on the grid that will lead to the
% correct predictive grid, i.e. compensates for forced grid movement
% INPUTS:
% measPdf - measurement pdf
% gridDimMeas - coordinates per dimension before rotation of measurement grid
% center - center of the measurement grid
% eigVect - eigvectors used to rotate the measurement grid
% Q - dynamics noise covariance matrix
% sFactor - wanted size of the grid based on the covariance of pdf
% nx - dimension of the state
% Npa - number of points per axis
% measMean - measurement pdf mean
% measVar - covariance of the measurement pdf
% OUTPUTS:
% measPdfNew - new interpolated measurement pdf
% gridDimNew - new coordinates per dimension before rotation
% gridStepNew - grid step of the new grid
% measGridNew - new measurement grid
% eigVectNew - eigenvectors used to rotate the new grid

%% Grid to interp on
mean = inv(F)*meanNew; %#ok<*MINV>
var = inv(F)*varNew*inv(F') + inv(F)*Q*inv(F');
[measGridNew, gridStepNew, gridDimNew, ~, eigVectNew] = gridCreation(mean,var,sFactor,nx,Npa); % create the grid to interp on

%% Interpolation
Fint = griddedInterpolant(gridDim,reshape(measPdf,Npa),"linear","none"); % interpolant
inerpOn = inv(eigVect)*(measGridNew - center); % Grid to inter on transformed to rotated space
measPdfNew = Fint(inerpOn')'; % Interpolation
measPdfNew(isnan(measPdfNew)) = 0; % Zeros for extrapolation, otherwise artifacts would appear

end

