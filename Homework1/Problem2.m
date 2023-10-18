%% Global variables
clear; clc; close all;

%% Physical parameters
% Density
rho_metal = 7000; % kg/m^3
rho_f = 1000; % fluid
rho = rho_metal - rho_f;

% Rod radius
r0 = 0.001;
% Young's modulus
Y = 1e9; % Using Y instead of E to avoid ambiguity
% Gravity
g = 9.8; % m/s^2
% Viscosity
visc = 1000; % Pa-s
% Total time
totalTime = 50; % seconds

% Rod length
RodLength = 0.1; % meter

% Time step size
dtt = [1e-4, 1e-3, 1e-2, 1e-1, 1e0]; % For Question 3
dt = dtt(3) % default = 1e-2

% Number of nodes
NN = [11, 21, 31, 41, 51]; % For Question 3
N = NN(2) % default = 21
ndof = N * 2; % # of DoF

% Discrete length
deltaL = RodLength / (N-1);

% Radius of spheres
R = zeros(N,1); % vector of size N
R(:) = deltaL/10;
midNode = (N+1)/2;
R(midNode) = 0.025;
 
% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - Bending stiffness
EA = Y * pi * r0^2; % Newton
% Geometry
nodes = zeros(N, 2); % Container of Geometry info

for c = 1:N
    nodes(c,1) = (c-1) * deltaL; % x coordinates of nodes
    % nodes(c,2) = 0; % y coordinates of nodes
end

% Mass matrix
M = zeros(ndof,ndof);
for k = 1:N
    M(2*k-1,2*k-1) = 4/3*pi*R(k)^3*rho_metal; % Mass for x_k
    M(2*k,2*k) = M(2*k-1,2*k-1); % Mass for y_k
end

% Viscous damping matrix
C = zeros(ndof,ndof);
for k=1:N
    C(2*k-1,2*k-1) = 6*pi*visc*R(k);
    C(2*k,2*k) = C(2*k-1,2*k-1);
end

% Gravity
W = zeros(ndof,1);
for k=1:N
    W(2*k) = -4/3*pi*R(k)^3*rho*g;
end

% Initial DOF vector
q0 = zeros(ndof,1); % initial position
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end

u0 = zeros(ndof,1); % initial velocity

% Tolerance
tol = EI / RodLength^2 * 1e-3;

% Number of time steps y-velocity of the middle node
Nsteps = round( totalTime / dt ) + 1;
all_mid_y = zeros( Nsteps, 1); % y-position of R2
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q0(N+1);
all_mid_v(1) = u0(N+1);

% Store information
Q = zeros(ndof, Nsteps);
Q(:,1) = q0;
U = zeros(ndof, Nsteps);
U(:,1) = u0;

% Time marching scheme

for c=2:Nsteps
 
%     fprintf('Time = %f\n', (c-1) * dt );
    
    q = q0; % Guess
    
    % Newton Raphson
    err = 10 * tol;
    while err > tol
        % Inertia
        f = M / dt * ( (q-q0) / dt - u0 );
        J = M / dt^2;
        %% Elastic forces
        % Linear spring between nodes k and k+1
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
        
        % Bending spring between k-1 and k+1. Rm: No bending at first or last node
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
        
        % Viscous force
        f = f + C * ( q - q0 ) / dt;
        J = J + C / dt;
        
        % Weight
        f = f - W;
        
        % Update
        q = q - J \ f;
        err = sum( abs(f) );
    end

    % Update
    u = (q - q0) / dt; % Velocity

    % New becomes Old
    q0 = q; Q(:,c) = q0;
    u0 = u; U(:,c) = u0;
    
%     % Plot to check the process
%     figure(1);
%     plot( q(1:2:end), q(2:2:end), 'ro-'); % plot 3 nodes
%     axis equal
%     xlabel('x (meter)')
%     ylabel('y (meter)')
%     drawnow
    
    % Store some info
    all_mid_y(c) = q(N+1);
    all_mid_v(c) = u(N+1);
end

%% Plotting
timeArray = (0:Nsteps-1) * dt;

% Q1:
figure(2); % Position and Velocity of Middle Node
subplot(1,2,1)
plot(timeArray, all_mid_y, 'r-','LineWidth',1.5);
xlabel('t (sec)');
ylabel('y (meter)');
xlim([0 50]);
% title('Position (along y-axis) of Middle Node');
grid on
subplot(1,2,2)
plot(timeArray, all_mid_v, 'b-','LineWidth',1.5);
xlabel('t (sec)');
ylabel('v [meter/sec]');
xlim([0 50]);
% title('Velocity (along y-axis) of Middle Node');
grid on

% Q2:
figure(3); % Shapes of structure
grid on;
hold on;
time_sep = [0, 0.05, 0.10, 1.0, 10.0, 50.0];
for i = 1:length(time_sep)
    ind = find(timeArray == time_sep(i));
    plot( Q(1:2:end,ind), Q(2:2:end,ind), 'o-','LineWidth',1.5); % plot N nodes
end
axis equal;
xlim([0 0.1]);
legend('t = 0','t = 0.05','t = 0.10','t = 1.0','t = 10.0','t = 50.0')
xlabel('x (meter)');
ylabel('y (meter)');
% title('Shapes of structure evolving');

% Q3:
figure(4); % Terminal Velocity vs. # of nodes and vs. time step size
subplot(1,2,1)
x1 = NN; % Data got from simulations
y1 = [-5.83707,-5.83427,-5.83375,-5.83357,-5.83348]*1e-3;
plot(x1, y1, 'ro-','LineWidth',1.5);
xlabel('Number of nodes');
ylabel('Terminal Velocity [meter/sec]');
ylim([-7 -4.5]*1e-3)
% title('Terminal velocity vs. the number of nodes');
grid on
subplot(1,2,2)
x2 = dtt; % Data got from simulations
y2 = [-5.83427,-5.83427,-5.83427,-5.83427,-5.83427]*1e-3;
semilogx(x2, y2, 'b*-','LineWidth',1.5);
xlabel('Time step size (sec)');
ylabel('Terminal Velocity [meter/sec]');
ylim([-7 -4.5]*1e-3)
% title('Terminal velocity vs. the time step size');
grid on

