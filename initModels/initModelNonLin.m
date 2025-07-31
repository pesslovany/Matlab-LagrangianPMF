function [model] = initModelNonLin(modelChoose)
%initModel initialize variables for estimation that are model dependent

switch modelChoose
    case 1 % NCV unknown turn rate
        model.kf = 100;
        model.nx = 5; % state dimension
        model.nz = 2; % measurement dimension
        model.dt = 1; % time step
        q = 0.01;
        model.Q = @(x) q*[(2*(x-sin(x)))./x.^3 (1-cos(x))./x.^2 zeros(1,1,size(x,3)) ((x-sin(x)))./x.^2 zeros(1,1,size(x,3));
            (1-cos(x))./x.^2 ones(1,1,size(x,3)) -((x-sin(x)))./x.^2 zeros(1,1,size(x,3)) zeros(1,1,size(x,3));
            zeros(1,1,size(x,3)) -((x-sin(x)))./x.^2 (2*(x-sin(x)))./x.^3 (1-cos(x))./x.^2 zeros(1,1,size(x,3));
            ((x-sin(x)))./x.^2 zeros(1,1,size(x,3)) (1-cos(x))./x.^2 ones(1,1,size(x,3)) zeros(1,1,size(x,3));
            zeros(1,1,size(x,3)) zeros(1,1,size(x,3)) zeros(1,1,size(x,3)) zeros(1,1,size(x,3)) repmat(deg2rad(0.01),1,1,size(x,3))]; % system noise
        model.u = zeros(model.nx,model.kf); % Input
        model.R = [deg2rad(0.01) 0 ;
            0 1]; % measurement noise covariance for both modes
        model.invR = inv(model.R);
        % PMF parameters
        model.Npa = [15 13 15 13 15];
        model.N = prod(model.Npa); % number of points - total
        model.sFactor = 4; % scaling factor (number of sigmas covered by the grid)
        model.meanV = [0;0]; % Mean values of components of meas noise
        model.wV = 1; % weights
        % Initial condition - Gaussian
        model.meanX0 = [100; 10; 100; 10; deg2rad(5)];% initial cond
        model.varX0 = diag([10 1 10 1 deg2rad(0.01)]); % initial cond variance
        model.f = @(x,u) [x(1,:) + (sin(x(5,:))./x(5,:)).*x(2,:) - ((1-cos(x(5,:)))./x(5,:)).*x(4,:);...
            cos(x(5,:)).*x(2,:) - sin(x(5,:)).*x(4,:);...
            ((1-cos(x(5,:)))/x(5,:)).*x(2,:) + x(3,:) + (sin(x(5,:))./x(5,:)).*x(4,:);
            sin(x(5,:)).*x(2,:) + cos(x(5,:)).*x(4,:);
            x(5,:)]; % model dynamics
        model.f_inv = @(y,u) [...
            y(1,:) - (sin(y(5,:))./y(5,:)).* (cos(y(5,:)).*y(2,:) + sin(y(5,:)).*y(4,:)) + ((1 - cos(y(5,:)))./y(5,:)).* (-sin(y(5,:)).*y(2,:) + cos(y(5,:)).*y(4,:)); ...
            cos(y(5,:)).*y(2,:) + sin(y(5,:)).*y(4,:); ...
            y(3,:) - ((1 - cos(y(5,:)))./y(5,:)).* (cos(y(5,:)).*y(2,:) + sin(y(5,:)).*y(4,:)) - (sin(y(5,:))./y(5,:)).* (-sin(y(5,:)).*y(2,:) + cos(y(5,:)).*y(4,:)); ...
            -sin(y(5,:)).*y(2,:) + cos(y(5,:)).*y(4,:); ...
            y(5,:) ...
            ];
        model.jacobian_f = @(x,u) [...
            1, (sin(x(5,:))./x(5,:)), 0, -(1 - cos(x(5,:)))./x(5,:), (x(2,:).*cos(x(5,:)) + x(4,:).*sin(x(5,:)))./x(5,:).^2 - (1 - cos(x(5,:)))./x(5,:).^2; ...
            0, cos(x(5,:)), 0, -sin(x(5,:)), -x(2,:).*sin(x(5,:)) - x(4,:).*cos(x(5,:)); ...
            0, (1 - cos(x(5,:)))./x(5,:), 1, (sin(x(5,:)))./x(5,:), -(x(2,:).*sin(x(5,:)))./x(5,:).^2 + (1 - cos(x(5,:)))./x(5,:).^2; ...
            0, sin(x(5,:)), 0, cos(x(5,:)), x(2,:).*cos(x(5,:)) - x(4,:).*sin(x(5,:)); ...
            0, 0, 0, 0, 1];
        model.x = mvnrnd(model.meanX0,model.varX0,1)';
        model.hfunct = @(x,v,k) [atan2(x(3,:), x(1,:)); sqrt(x(1,:).^2 + x(3,:).^2)] + v;
        for k = 1:model.kf-1
            w = mvnrnd(zeros(model.nx,1),model.Q(model.x(5,k)), 1)';
            model.x(:,k+1) = model.f(model.x(:,k),model.u(:,k)) + w;
        end
        % measurement generation
        model.z = model.hfunct(model.x,0,0)+sqrt(model.R)*randn(model.nz,width(model.x));
        model.V.pdf = gmdistribution(model.meanV',repmat(model.R,1,1,1),model.wV); % Meas noise pdf
    case 2 % Nonholonomic robot with bearing‐only measurements
        landmarks = [  0, 10,  0, 10;
            0,  0, 10, 10 ];
        model.kf      = 200;
        model.nx      = 3;                                    % [x; y; θ]
        model.nz      = size(landmarks,2);                    % one bearing per landmark
        model.dt      = 1;                                    % time step
        q             = 1e-2;                                 % much lower process noise
        model.Q       = @(x) q * eye(model.nx);               % simple isotropic process noise
        model.u       = repmat([1.0; 0.05; 0], 1, model.kf);  % forward speed = 1 m/s, turn rate = 0.05 rad/s
        model.R       = (deg2rad(4)^2) * eye(model.nz);     % slightly less measurement noise
        model.invR    = inv(model.R);
        model.Npa     = [71 71 71];                           % PMF grid sizes
        model.N       = prod(model.Npa);
        model.sFactor = 4;
        model.meanV   = zeros(model.nz,1);
        model.wV      = 1;
        model.meanX0  = [0; 0; deg2rad(1)];                   % [x0; y0; θ0]
        model.varX0   = diag([1 1 deg2rad(0.005)]);           % tighter initial covariance
        model.f       = @(x,u) [ ...
            x(1,:) + model.dt * u(1,:) .* cos(x(3,:)); ...
            x(2,:) + model.dt * u(1,:) .* sin(x(3,:)); ...
            x(3,:) + model.dt * u(2,:) ];
        model.f_inv   = @(y,u) [ ...
            y(1,:) - model.dt * u(1,:) .* cos(y(3,:)); ...
            y(2,:) - model.dt * u(1,:) .* sin(y(3,:)); ...
            y(3,:) - model.dt * u(2,:) ];
        model.jacobian_f = @(x,u) [ ...
            1, 0, -model.dt * u(1,:) .* sin(x(3,:)); ...
            0, 1,  model.dt * u(1,:) .* cos(x(3,:)); ...
            0, 0,  1 ];
        model.x       = mvnrnd(model.meanX0, model.varX0, 1)'; % propagate initial state
        model.hfunct  = @(x,v,k) wrapToPi( ...
            atan2( ...
            repmat(landmarks(2,:)',1,size(x,2)) - repmat(x(2,:),size(landmarks,2),1), ...
            repmat(landmarks(1,:)',1,size(x,2)) - repmat(x(1,:),size(landmarks,2),1) ...
            ) ...
            - repmat(x(3,:),size(landmarks,2),1) ...
            ) + v;

        for k = 1:model.kf-1
            w = mvnrnd( zeros(model.nx,1), model.Q(model.x(3,k)), 1 )';
            model.x(:,k+1) = model.f(model.x(:,k),model.u(:,k)) + w;
        end
        model.z       = model.hfunct(model.x, 0, 0) + sqrt(model.R) * randn(model.nz, size(model.x,2));
        model.V.pdf   = gmdistribution( model.meanV', repmat(model.R,1,1,1), model.wV );
end

end