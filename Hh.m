function J = Hh(x, vysky)
    % Compute the gradient for the Jacobian matrix using central differences
    dx = (eps(x(:, 1))).^(1/3); % Small step sizes for finite difference
    dx1 = [dx(1); 0];           % Step in the x direction
    dx2 = [0; dx(2)];           % Step in the y direction

    % Compute function values for central difference directly using vysky
    J11 = vysky(x(1,:) + dx1(1), x(2,:)); % f(x + dx1)
    J12 = vysky(x(1,:) - dx1(1), x(2,:)); % f(x - dx1)
    J21 = vysky(x(1,:), x(2,:) + dx2(2)); % f(x + dx2)
    J22 = vysky(x(1,:), x(2,:) - dx2(2)); % f(x - dx2)

    % Central difference for gradients
    J1 = (J11 - J12) / (2 * dx1(1)); % Gradient with respect to x
    J2 = (J21 - J22) / (2 * dx2(2)); % Gradient with respect to y

    J3 = ones(size(J1));

    % Reshape the Jacobian into the desired format
    J = reshape([J1; J2; J3], 1, 3, []);

end
