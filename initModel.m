function [model] = initModel(modelChoose, souradniceGNSS, hB, vysky)
%initModel initialize variables for estimation that are model dependent

% Time steps
timeSteps = 1:1:length(souradniceGNSS);
% Last time
endTime = length(timeSteps);


switch modelChoose
    case 3 % 3D
        % System parameters
        model.nx = 3; % state dimension
        model.nz = 1; % measurement dimension
        model.dt = 1; % time step
        model.Q = blkdiag(36*eye(2),0.05); % system noise
        model.u = [ [diff(souradniceGNSS(1:2,timeSteps),1,2) [0;0]] + sqrt(model.Q(1:2,1:2))*randn(2,endTime); zeros(1,endTime)] ; % Input
        model.invQ = inv(model.Q);
        model.R = 3; % measurement noise covariance for both modes
        model.invR = inv(model.R);
        % PMF parameters
        model.Npa = [71 71 21]; % number of points per axis
        model.N = prod(model.Npa); % number of points - total
        model.noPart = 200000; % number of particles for PF
        model.sFactor = 5; % scaling factor (number of sigmas covered by the grid)
        model.meanV = 0; % Mean values of components of meas noise
        model.wV = 1; % weights
        model.x = [souradniceGNSS(1:2,timeSteps);hB-vysky(souradniceGNSS(1,timeSteps),souradniceGNSS(2,timeSteps))];
        % Initial condition - Gaussian
        model.meanX0 = [souradniceGNSS(1:2,timeSteps(1)); 0];% initial cond
        model.varX0 = [25 0 0;
            0 25 0;
            0 0 10]; % initial cond variance
        model.F = [1 0 0;
            0 1  0;
            0 0 0.96925]; % model dynamics
        model.hfunct = @(x,v,k) [vysky(x(1,:),x(2,:))] + x(3,:) + v; % measurement equation
        % measurement generation
        model.z = model.hfunct(model.x,0,0)+sqrt(model.R)*randn(model.nz,length(model.x));
        model.z(1,:) = hB; % add real measurements
        model.V.pdf = gmdistribution(model.meanV',repmat(model.R,1,1,1),model.wV); % Meas noise pdf
        % UKF Params
        model.ukfOn = 1; % 0 - disabled, 1 - enabled
        model.kappa = 1;
        model.SPnum = 2*model.nx+1;
    case 4 % 4D
        model.nx = 4; % state dimension
        model.nz = 3; % measurement dimension
        model.u = zeros(model.nx,endTime);
        model.dt = 1; % time step
        model.q = 1; % noise parameter
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
        model.noPart = floor(1.8*model.N); % number of particles for PF
        model.sFactor = 6; % scaling factor (number of sigmas covered by the grid)
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
        model.ukfOn = 1; % 0 - disabled, 1 - enabled
        model.kappa = 1;
        model.SPnum = 2*model.nx+1;
    otherwise
        disp('Error: unsupported model')
end

end