% laplace_relax.m
%-------------------------------------------------------------------------------
% Solve the Laplace equation with Dirichlet BCs using relaxation
% (Jacobi or Gauss-Seidel)
%-------------------------------------------------------------------------------
% Clear memory and show only a few digits
clear('all');
format('short');

%-------------------------------------------------------------------------------
% Select which method to use
%-------------------------------------------------------------------------------
whatMethod = 'Jacobi';
% whatMethod = 'GaussSeidel';
% whatMethod = 'SOR';
fprintf(1,'Using %s Relaxation\n',whatMethod);

%-------------------------------------------------------------------------------
% Set parameters
%-------------------------------------------------------------------------------
h = 0.05;         % Spatial step
maxIterations = 1e+4;  % Maximum number of iterations
minDiff = 1e-6;  % Convergence criterion
frameSkip = 1;  % Only plot every frameSkip iterations
frameUpdateLag = 0; % Slow down animation

% Set up vectors of discretized x and y values
x = 0:h:1;
y = x;
L = length(x);

%-------------------------------------------------------------------------------
% Define initial condition (phi)
phi_new = zeros(L);
% Impose BCs
phi_new(2:L-1,L) = 1;
phi = phi_new;

%-------------------------------------------------------------------------------
% Set over-relaxation factor, omega for SOR
if strcmp(whatMethod,'SOR')
    %omega=1.5;
    %omega=1.2;
    %omega=1;
    % Analytic result for optimal value
    omega_opt = 2/(1+sin(pi/L));
    fprintf(1,'Optimal omega: %g\n',omega_opt);
    omega = omega_opt;
end

%-------------------------------------------------------------------------------
% Set up figure for plotting
f = figure(1);
f.Color = 'w';
set_color_map()
h_surface = surfc(x,y,phi');
xlabel('x');
ylabel('y');
zlabel('Potential \phi');
shading('faceted');

%-------------------------------------------------------------------------------
% Relaxation loop
for iter = 1:maxIterations

    switch whatMethod
    case 'Jacobi'
        % Loop over interior points and update.
        % Note the boundary values are preserved.
        for j = 2:L-1
            for l = 2:L-1
                phi_new(j,l) = 0.25*(phi(j-1,l) + phi(j+1,l) + ...
                            phi(j,l-1) + phi(j,l+1));
            end
        end

    case 'GaussSeidel'
        % Keep this for the convergence test
        old_phi = phi;

        % Loop over interior points: update in place, and preserve BCs
        for j = 2:L-1
            for l = 2:L-1
                phi(j,l) = 0.25*(phi(j-1,l) + phi(j+1,l) + ...
                                phi(j,l-1) + phi(j,l+1));
            end
        end

        % Update:
        phi_new = phi;
        phi = old_phi;

    case 'SOR'
        % Successive overrelaxation (SOR)
        % Keep this for the convergence test
        old_phi = phi;

        % Loop over interior points: update in place, and preserve BCs
        for j = 2:L-1
            for l = 2:L-1
                gs = 0.25*(phi(j-1,l) + phi(j+1,l) + ...
                                phi(j,l-1) + phi(j,l+1));
                phi(j,l) = (1-omega)*phi(j,l) + omega*gs;
            end
        end

        % Update:
        phi_new = phi;
        phi = old_phi;
    end

    % Plot updated solution
    if rem(iter,frameSkip) == 0
        h_surface(1).ZData = phi_new';
        title(sprintf('Iteration: %u',iter));
        drawnow();
        pause(frameUpdateLag)
    end

    % Break if change is suitably small
    diff = max(abs(phi_new - phi));
    if diff < minDiff
        disp('Relaxation converged...');
        break;
    end

    % To be updated in the next iteration:
    phi = phi_new;
end

% Display number of iterations
fprintf(1,'%u Iterations\n',iter);

%-------------------------------------------------------------------------------
% Determine the electric field by differencing. Note Matlab's
% identification of the first index in a matrix with
% y, wherever a functional dependence z = z(x,y) is implied.
[Ey, Ex] = gradient(phi,h);
Ex = -Ex;
Ey = -Ey;

%-------------------------------------------------------------------------------
% Also plot the electric field
f2 = figure(2);
f2.Color = 'w';
hold('on');
contourf(x,y,phi');
scale = 10;
quiver(x,y,Ex',Ey',scale,'w');
xlabel('x');
ylabel('y');
title('Electric field {\bf E} and contours of \phi');
axis([0 1 0 1]);
hold('off');
set_color_map()
