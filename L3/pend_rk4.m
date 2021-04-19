% pend_rk4.m
% Solve the nonlinear pendulum problem using RK4

% Clear memory and show only a few digits
clear('all');
format('short');

%-------------------------------------------------------------------------------
% Time grid parameters:
tau = 0.025; % Time step
T = 10; % Total integration time
numSteps = ceil(T/tau); % Number of time steps

%-------------------------------------------------------------------------------
% Initial conditions
theta1 = 180; % Initial angle in degrees
theta = theta1*pi/180;
x = [theta 0];

%-------------------------------------------------------------------------------
% Preallocate:
time = zeros(numSteps,1);

%-------------------------------------------------------------------------------
% Fourth-order Runge-Kutta integration
figure(1);
for n = 1:numSteps
    % Time
    time(n) = (n-1)*tau;

    % One step of RK4
    f1 = rhs_pend(x);
    f2 = rhs_pend(x+0.5*tau*f1);
    f3 = rhs_pend(x+0.5*tau*f2);
    f4 = rhs_pend(x+tau*f3);
    x = x+tau*(f1+2*f2+2*f3+f4)/6;

    % Co-ordinates of the pendulum bar
    xPend = [0  sin(x(1))];
    yPend = [0 -cos(x(1))];

    % Animate the pendulum motion
    plot(xPend,yPend,'o-')
    title(sprintf('Time: %5.3f',time(n)+tau))
    axis('equal')
    axis([-1 1 -1 1])
    pause(0.05);
    drawnow;
end
