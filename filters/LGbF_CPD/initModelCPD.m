function model = initModelCPD(souradniceGNSS, hB, vysky)

% 4D NCV with measurement altituted and velocity in body frame
% State and measurement dimensions
model.nx = 4;
model.nz = 3;

% Time and input
model.timeSteps = 1:size(souradniceGNSS, 2);
model.endTime = numel(model.timeSteps);
model.dt = 1;
model.u = zeros(model.nx, model.endTime);

% System noise
q = 4;
dt = model.dt;
model.Q = q * [dt^3/3 0 dt^2/2 0;
    0 dt^3/3 0 dt^2/2;
    dt^2/2 0 dt 0;
    0 dt^2/2 0 dt];
model.invQ = inv(model.Q);

% Measurement noise
model.R = diag([3, 1, 1].^2);
model.invR = inv(model.R);

% PMF grid parameters
model.Npa = [51 51 41 41];
model.sFactor = 6;
model.sFactorCPD = 6;

% Measurement noise PDF
meanV = [0; 0; 0];
wV = 1;
model.V.pdf = gmdistribution(meanV', repmat(model.R, 1, 1, 1), wV);
model.Vsplit(1).pdf = gmdistribution(meanV(1), model.R(1,1), 1);
model.Vsplit(2).pdf = gmdistribution(meanV(2), model.R(2,2), 1);
model.Vsplit(3).pdf = gmdistribution(meanV(3), model.R(3,3), 1);

% Initial condition
velikost = [diff(souradniceGNSS(1:2, model.timeSteps), 1, 2), [0; 0]];
model.meanX0 = [souradniceGNSS(1:2, model.timeSteps(1)); velikost(:, 1)];
model.varX0 = [25 0 0 0;
    0 25 0 0;
    0 0 0.5 0;
    0 0 0 0.5];

% Dynamics (NCV)
model.F = [1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1];

% Measurement equation
model.hfunct = @(x, v, k) [
    vysky(x(1,:), x(2,:));
    [ (x(3,:)./sqrt(x(3,:).^2 + x(4,:).^2)).*x(3,:) - (x(4,:)./sqrt(x(3,:).^2 + x(4,:).^2)).*x(4,:);
    (x(4,:)./sqrt(x(3,:).^2 + x(4,:).^2)).*x(3,:) + (x(3,:)./sqrt(x(3,:).^2 + x(4,:).^2)).*x(4,:)
    ]
    ] + v;

model.hFunctSplit{1} = @(x) vysky(x(1,:), x(2,:));
model.hFunctSplit{2} = @(x) (x(1,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(1,:) - (x(2,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(2,:);
model.hFunctSplit{3} = @(x) (x(2,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(1,:) + (x(1,:)./sqrt(x(1,:).^2 + x(2,:).^2)).*x(2,:);

model.mapping.Combs = {[1 2], [3 4]};
model.mapping.indexed = [1 2 2];

% True state and measurements
model.x = [souradniceGNSS(1:2, model.timeSteps); velikost];
model.z = model.hfunct(model.x, 0, 0);
model.z(1,:) = hB;
model.z(2:3,:) = model.z(2:3,:) + chol(model.R(2:3,2:3)) * randn(size(model.z(2:3,:)));

model.NpaCPD = model.Npa;
model.N = prod(model.Npa);


end