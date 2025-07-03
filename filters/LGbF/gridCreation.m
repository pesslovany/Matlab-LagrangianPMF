function [gridOut, gridStep, gridCoordPerDim, meanIn, eigVect, gridBound] = gridCreation(meanIn,covIn,sFactor,nx,Npa)
%gridCreation creates a grid based on first two moments
% INPUTS:
% meanIn - center of the grid
% covIn - covariance giving the size and rotation of the grid
% sFactor - size of the grid based on covariance
% nx - dimension of the grid
% Npa - vector, number of points per axe
% OUTPUTS:
% gridOut - created grid
% gridStep - grid step size per dim
% gridCoordPerDim - grid coordinates per dimension before rotation
% meanIn - center of the created grid
% eigVect - eigenvectors that were used to rotate the grid

[eigVect,eigVal] = eig(covIn); % eigenvalue and eigenvectors, for setting up the grid
eigVal = diag(eigVal);
gridBound = sqrt(eigVal)*sFactor; % Boundaries of grid

% Ensure the grid steps are in the right order
[~,I] = sort(diag(covIn));
[~,I] = sort(I);

[pom,Ipom] = sort(gridBound);
gridBound = pom(I);

pom2 = eigVect(:,Ipom);
eigVect = pom2(:,I);

%Creation of grid
gridCoordPerDim = cell(nx,1);
gridStep = nan(nx,1);
for ind3 = 1:1:nx 
    gridCoordPerDim{ind3,1} = linspace(-gridBound(ind3),gridBound(ind3),Npa(ind3)); %New grid with middle in 0 in one coordinate
    gridStep(ind3,1) = abs(gridCoordPerDim{ind3,1}(1)-gridCoordPerDim{ind3,1}(2)); %Grid step
end
gridOut = eigVect*combvec(gridCoordPerDim) + meanIn; %Grid rotation by eigenvectors and traslation to the mean

end

