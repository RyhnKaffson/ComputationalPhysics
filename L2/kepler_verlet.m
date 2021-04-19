% kepler_verlet.m
% Motion under a central force, using the Verlet method

% Clear memory and only show a few digits
clear('all');
format('short');

% Time step (non-dim.)
tau = 0.05;

% Initial position (non-dim.) - this should be fixed
pos = [1 0];

% Initial velocity (non-dim.) - vary the y-component
vel = [0 0.3];

% Total integration time
T = 810;

% Number of integration steps
numSteps = ceil(T/tau);

% Avoid plotting every step
skip = ceil(numSteps/50);

% Calculate trajectory from analytic solution. See Appendix
% of Week 2 lecture 1.
[xan,yan] = kepler_analytic(vel,T);

% Preallocate vectors for speed:
x = zeros(numSteps,1);
y = zeros(numSteps,1);
time = zeros(numSteps,1);
energy = zeros(numSteps,1);

%-------------------------------------------------------------------------------
% Verlet method integration
%-------------------------------------------------------------------------------
figure(1);
for n = 1:numSteps
    % Store position and time for plotting: time-step n
    x(n) = pos(1);
    y(n) = pos(2);
    time(n) = (n-1)*tau;

    % Plot numerical and analytic solution
    if rem(n,skip)==0
        plot(x,y,'g-',pos(1),pos(2),'ko',xan,yan,'b',0,0,'ro')
        title(['Time: ',num2str(time(n))]); xlabel('x'); ylabel('y');
        axis equal; % Preserve aspect ratio
        drawnow(); % Draw immediately
    end

    % Calculate radial position and acceleration at step n
    r = norm(pos);
    accel = -pos/r^3;

    % Verlet Method: position - prev is step n-1, next is step n+1
    if n==1
        % Get started with a midpoint method step
        next = pos+tau*vel+0.5*tau^2*accel;
    else
        next = 2*pos-prev+tau^2*accel;
    end

    % Verlet method: velocity at step n
    if n > 1
        vel = (next-prev)/(2*tau);
    end
    speed = norm(vel);

    % Calculate total energy at step n and store
    energy(n) = 0.5*speed^2-1/r;

    % Update
    prev = pos;
    pos = next;
end

% Plot energy versus time
figure(2);
plot(time,energy);
xlabel('Time (non-dim.)');
ylabel('Total energy (non-dim.)');
