function [model] = initModelCPD(souradniceGNSS, hB, vysky)
%initModel initialize variables for estimation that are model dependent

% Time steps
timeSteps = 1:1:length(souradniceGNSS);
% Last time
endTime = length(timeSteps);



model.nx = 4; % state dimension
model.nz = 3; % measurement dimension
model.u = zeros(model.nx,endTime);
model.dt = 1; % time step
model.q = 4; % noise parameter
model.Q = model.q*[ model.dt^3/3 0 model.dt^2/2 0;
    0 model.dt^3/3 0 model.dt^2/2;
    model.dt^2/2 0 model.dt 0;
    0 model.dt^2/2 0 model.dt]; % system noise
model.invQ = inv(model.Q);
model.R = diag([3, 1, 1].^2); % measurement noise covariance for both modes %2
model.invR = inv(model.R);
% PMF parameters
model.Npa = [51 51 41 41]; % number of points per axis
model.N = prod(model.Npa); % number of points - total
model.noPart = model.N;%floor(1.8*model.N); % number of particles for PF
model.sFactor = 5; % scaling factor (number of sigmas covered by the grid)
model.essThrd = (2/3)*model.noPart; % effective sample size threshold for PF
model.meanV = [0; 0; 0]; % Mean values of components of meas noise
model.wV = 1; % weights
model.V.pdf = gmdistribution(model.meanV',repmat(model.R,1,1,1),model.wV); % Meas noise pdf
velikost = [diff(souradniceGNSS(1:2,timeSteps),1,2) [0;0]];
% Initial condition - Gaussian
model.meanX0 = [souradniceGNSS(1:2,timeSteps(1));velikost(:,1)];% initial cond
model.varX0 = [25 0 0 0;
    0 25 0 0 ;
    0 0 0.5 0;
    0 0 0 0.5]; % initial cond variance
% Dynamics - known turn rate
model.F = [1 0 model.dt 0;
    0 1 0 model.dt;
    0 0 1 0;
    0 0 0 1];
model.hfunct = @(x,v,k) [vysky(x(1,:),x(2,:)); [(x(3,:)./sqrt(x(3,:).^2+x(4,:).^2)).*x(3,:)-(x(4,:)./sqrt(x(3,:).^2+x(4,:).^2)).*x(4,:); (x(4,:)./sqrt(x(3,:).^2+x(4,:).^2)).*x(3,:)+(x(3,:)./sqrt(x(3,:).^2+x(4,:).^2)).*x(4,:)]  ] + v; % measurement equation
model.x = [souradniceGNSS(1:2,timeSteps); velikost]; % state
model.z = model.hfunct(model.x,0,0); % Measurements
model.z(1,:) = hB; % add real data as measurement
model.z(2:3,:) = model.z(2:3,:)+chol(model.R(2:3,2:3))*randn(size(model.z(2:3,:)));
% UKF Params
model.kappa = 1;
model.SPnum = 2*model.nx+1;
% RPBF Params
model.xnid = [1, 2]; % indices for nonlinear states (RBPF)
model.xlid = [3, 4]; % indices for linear states (RBPF)
model.nxn = length(model.xnid);
model.nxl = length(model.xlid);
model.hJacobLfunct = @(x) [0, 0;
    (x(3)^3 + 3*x(3)*x(4)^2)/sqrt((x(3)^2 + x(4)^2)^3), -(3*x(3)^2*x(4) + x(4)^3)/sqrt((x(3)^2 + x(4)^2)^3);
    2*x(4)^3/sqrt((x(3)^2 + x(4)^2)^3),                 2*x(3)^3/sqrt((x(3)^2 + x(4)^2)^3)]; % Jacobian matrix of the observation model for the linear state (RBPF)
model.noPartRbpf = 20000;    % number of particles for RBPF
model.essThrdRbpf = (2/3)*model.noPart; % effective sample size threshold for RBPF

% CPD
model.V.pdf = gmdistribution(model.meanV', repmat(model.R, 1, 1, 1), model.wV);
model.Vsplit(1).pdf = gmdistribution(model.meanV(1), model.R(1,1), 1);
model.Vsplit(2).pdf = gmdistribution(model.meanV(2), model.R(2,2), 1);
model.Vsplit(3).pdf = gmdistribution(model.meanV(3), model.R(3,3), 1);

model.hFunctSplit{1} = @(x) vysky(x(1,:), x(2,:));
model.hFunctSplit{2} = @(x) (x(1,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(1,:) - (x(2,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(2,:);
model.hFunctSplit{3} = @(x) (x(2,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(1,:) + (x(1,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(2,:);

model.mapping.Combs = {[1 2], [3 4]};
model.mapping.indexed = [1 2 2];

model.NpaCPD = [101 101 101 101];
% Max rank of CPD
model.Rdef = 10;
model.Q = diag(diag(model.Q));
model.sFactorCPD = 6;


end
