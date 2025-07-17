function [measPdf, GridDelta, xOld, predGridAdvect, gridDimOld, ...
          Ppold, gridBound, measGridNew, predVarPom] = ...
          measPdfInterpSpect(measPdf, GridDelta, measMean, measVar, ...
                          k, u, F, Q, sFactor, nx, Npa, ...
                          gridDim, gridRotation, gridCenter)
%GBFADVECTINTERP

predVarPom       = F*measVar(:,:,k)*F' + Q;
xOld            = F*measMean(:,k) + u(:,k);

[predGridAdvect, GridDelta(:,k+1), gridDimOld, ~, Ppold, gridBound] = ...
    gridCreation(xOld, predVarPom, sFactor, nx, Npa);

measGridNew     = F \ (predGridAdvect - u(:,k));
GridDelta(:,k)  = F \ GridDelta(:,k+1);

Fint             = griddedInterpolant(gridDim, reshape(measPdf, Npa), "linear", "none");
inerpOn          = inv(gridRotation) * (measGridNew - gridCenter);
measPdf         = Fint(inerpOn');
measPdf(isnan(measPdf)) = 0;

end