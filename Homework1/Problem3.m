%% Global variables
clear; clc; close all;

%% Physical parameters
% Density
rho_metal = 2700; % kg/m^3

% Rod radius
r = 0.011; % inner radius
R = 0.013; % outer radius

% Young's modulus
Y = 70*1e9; % Using Y instead of E to avoid ambiguity
% Gravity
g = 9.8; % m/s^2
% Total time
totalTime = 1; % seconds

% Rod length
RodLength = 1; % meter

% Time step size
dt = 1e-2;

% Number of nodes
N = 50;
ndof = N * 2; % # of DoF

% Discrete length
deltaL = RodLength / (N-1);

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * (R^4-r^4) / 4; % Nm^2 - Bending stiffness
EA = Y * pi * (R^2-r^2); % Newton

% Mass matrix
dm = pi * (R^2-r^2) * RodLength * rho_metal / (N - 1);
M = eye(ndof) * dm;

% Tolerance
tol = EI / RodLength^2 * 1e-3;

% Number of time steps y-velocity of the middle node
Nsteps = round( totalTime / dt ) + 1;

% Fixed and Free DoFs
fixedDOF = [1; 2; ndof];
freeDOF = 3:ndof-1;
boundaryConditionVector = [0; 0; 0];

% "Gravity" - P External Force
PP0 = 2e3:1e3:2e4;
YY_max = zeros(1,length(PP0));
yy_max = zeros(1,length(PP0));

for i = 1:1 % 1:length(PP0) for Q2
% Geometry
nodes = zeros(N, 2); % Container of Geometry info
for c = 1:N
    nodes(c,1) = (c-1) * deltaL; % x coordinates of nodes
    % nodes(c,2) = 0; % y coordinates of nodes
end

% Find the closest node
Point = 0.75; % length from the left side
[~,ind_P] = min(abs(nodes(:,1) - Point));
P0 = PP0(i); %%%%%%%%%%%%% default = -2000
P = zeros(ndof,1);
P(2*ind_P) = -P0;

% Initial DOF vector
q0 = zeros(ndof,1); % initial position
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end

u0 = zeros(ndof,1); % initial velocity

% Store information
Q = zeros(ndof, Nsteps);
Q(:,1) = q0;
Y_max = zeros(1, Nsteps);

% Time marching scheme

for c=2:Nsteps
 
%     fprintf('Time = %f\n', (c-1) * dt );
    
    q = q0; % Guess
    q(fixedDOF) = boundaryConditionVector;
    
    % Newton Raphson
    err = 10 * tol;
    while err > tol

        q_free = q(freeDOF);

        % Inertia
        f = M / dt * ( (q-q0) / dt - u0 );
        J = M / dt^2;
        %% Elastic forces
        % Stretching energy between nodes k and k+1
        for k = 1:N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            ind = 2*k-1:2*k+2;% Include all related nodes' index
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end
        
        % Bending energy between k-1 and k+1. Rm: No bending at first or last node
        for k = 2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            ind = 2*k-3:2*k+2;
            f(ind) = f(ind) + dF;
            J(ind,ind) = J(ind,ind) + dJ;
        end
        
        % Weight
        f = f - P;
        
        %%%%%%%%%%%%
        f_free = f(freeDOF);
        J_free = J(freeDOF,freeDOF);
        
        % Update
        dq_free = J_free \ f_free;
        q_free = q_free - dq_free;

        err = sum( abs(f_free) ); % only focus on free parts

        q(freeDOF) = q_free;

    end

    % Update
    u = (q - q0) / dt; % Velocity

    % New becomes Old
    q0 = q; Q(:,c) = q0;
    u0 = u;

    Y_max(c) = min(q0);
    
%     % Plot to check the process
%     figure(1);
%     plot( q(1:2:end), q(2:2:end), 'ro-');
%     axis equal
%     xlabel('x (meter)')
%     ylabel('y (meter)')
%     drawnow

end

%% Plotting
timeArray = (0:Nsteps-1) * dt;

% Q1:
figure(2); % Maximum vertical displacement
plot(timeArray, Y_max, 'r-','LineWidth',1.5);
xlabel('t (sec)');
ylabel('y_{max} (meter)');
xlim([0 1]);
% title('Maximum vertical displacement of the beam');
grid on

% Euler Beam Theory:
temp = min(Point, RodLength - Point);
y_max = -P0*temp*(RodLength^2-temp^2)^1.5/(9*sqrt(3)*EI*RodLength);

YY_max(i) = Y_max(c);
yy_max(i) = y_max;

end

% Q2:
figure(3); % P vs. y_max
grid on;
hold on;
plot( PP0, YY_max, 'ro-','LineWidth',1.5); % plot Simulation results
plot( PP0, yy_max, 'b*-','LineWidth',1.5); % plot Theory results
legend('Simulation Results','Theory Results')
xlabel('P (Newton)');
ylabel('y_{max} (meter)');
% title('Simulation and theory difference');
