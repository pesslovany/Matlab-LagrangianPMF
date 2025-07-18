function [measPdfNew,gridDimNew, gridStepNew, measGridNew, eigVectNew, measAdvWant] = measPdfInterpNonLin(measPdf,gridDim,center,eigVect,Q,sFactor,nx,Npa,measMean,measVar,f,invf,jacobianf,u)
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
measAdvWant = f(measMean,u);
varAdvWant = jacobianf(measMean,u)*measVar*jacobianf(measMean,u)' + Q(measMean(3));

% Advected grid
[predGrid, predGridDelta, gridDimNew, ~, eigVectNew] = gridCreation(measAdvWant,varAdvWant,sFactor,nx,Npa); % create the grid to interp on
measGridNew = invf(predGrid,u); % Interpolating grid

filtGridDelta_nonlin = abs(invf(predGrid+predGridDelta,u)-measGridNew);
gridStepNew = prod(filtGridDelta_nonlin)';

%% Interpolation
Fint = griddedInterpolant(gridDim,reshape(measPdf,Npa),"linear","none"); % interpolant
inerpOn = inv(eigVect)*(measGridNew - center); % Grid to inter on transformed to rotated space
measPdfNew = Fint(inerpOn')'; % Interpolation
measPdfNew(isnan(measPdfNew)) = 0; % Zeros for extrapolation, otherwise artifacts would appear

measPdfNew = measPdfNew./sum(measPdfNew.*gridStepNew')'; % Normalizaton (theoretically not needed)


end

