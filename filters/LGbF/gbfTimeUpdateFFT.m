function [predPdf,predGrid,gridStep] = gbfTimeUpdateFFT(F,measPdf,measGrid,gridStep,k,Npa,invQ,predDenDenomW,nx,u)
%PMFUPDATEFFT time update in the Lagrangian formulation
%INPUTS:
% F - model dynamics matrix
% measPdf - measurement pdf
% measGrid - measurement grid
% gridStep - measurement grid step per dimension (and all old grid steps)
% k - time step
% Npa - number of points per axis
% invQ - inversion of dynamics noise
% predDenDenomW - constant of the dynamics gaussian noise function 
% nx - dimension of the state
% u - input to the model dynamics
%OUTPUTS:
% predPdf - predictive pdf
% predGrid - predictive grid
% gridStep - all grids steps with added predictive grid step


% Pred Grid
predGrid = F*measGrid + u; % Predictive grid
gridStep(:,k+1) = F*gridStep(:,k); % Predictive grid step size

measPdfDotDeltas = (measPdf*prod(gridStep(:,k))); % measurement PDF * measurement PDF step size
mesPdfDotDeltasTens = reshape(measPdfDotDeltas,Npa); % Into tensor space

halfGridInd = ceil(length(predGrid)/2); % Index of middle point of predictive grid
distPoints = (predGrid(:,halfGridInd)'-(predGrid)'); % Distance of transformed grid points to the new grid points
convKer = ((exp(sum(-0.5*distPoints*invQ.*distPoints,2)))/predDenDenomW)';% Convolution kernel values
convKerTens = reshape(convKer,Npa); % Convolution kernel values in tensor format

dims = 1:1:nx;
% Will be used to truncate the padding need the do the convolution
ifun = @(m,n) ceil((n-1)/2)+(1:m);
subs(1:ndims(mesPdfDotDeltasTens)) = {':'};

for dim=dims % FFT over all dimensions
    % compute the FFT length with padding
    l = Npa(dim)+Npa(dim)-1;
    mesPdfDotDeltasTens = fft(mesPdfDotDeltasTens,l,dim); % FFT measurement pdf
    convKerTens = fft(convKerTens,l,dim); % FFT of the kernel
    subs{dim} = ifun(Npa(dim),Npa(dim)); % Padding indices
end

% Use this in case the mex file not working
% mesPdfDotDeltasTens = mesPdfDotDeltasTens.*convKerTens;

% Convolution in frequency domain (inplaceprod saves memory)
inplaceprod(mesPdfDotDeltasTens, convKerTens);

% IFFT back to space domain
for dim=dims
    mesPdfDotDeltasTens = ifft(mesPdfDotDeltasTens,[],dim);
end
% Make sure the result is real and delete the padding
predDensityProb2cub = real(mesPdfDotDeltasTens(subs{:}));

predPdf = reshape(predDensityProb2cub,length(predGrid),1); % back to vector for easier manipulation
predPdf = predPdf./(sum(predPdf)*prod(gridStep(:,k+1)))'; % Normalizaton (theoretically not needed)
predPdf(predPdf<0) = 0; % FFT approach sometimes can produce very small negative values

end

