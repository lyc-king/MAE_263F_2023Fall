%% Global variables
clear; clc; close all;

%% Physical parameters
% Number of nodes
N = 3;
ndof = N * 2; % # of DoF

% Rod length
RodLength = 0.1; % meter

% Discrete length
deltaL = RodLength / (N-1);

% Radius of spheres (m)
R1 = 0.005; % 0.025 for Q3
R2 = 0.025;
R3 = 0.005; % 0.025 for Q3

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
totalTime = 10; % seconds

% Utility quantities
ne = N - 1; % Number of edges
EI = Y * pi * r0^4 / 4; % Nm^2 - Bending stiffness
EA = Y * pi * r0^2; % Newton
% Geometry
nodes = zeros(N, 2); % Container of Geometry info

for c = 1:N
    nodes(c,1) = (c-1) * deltaL; % x coordinates of nodes at beginning
%     nodes(c,2) = 0; % y coordinates of nodes
end

% Mass matrix
M = zeros(ndof,ndof);
M(1,1) = 4/3*pi*R1^3*rho_metal;
M(2,2) = M(1,1);
M(3,3) = 4/3*pi*R2^3*rho_metal;
M(4,4) = M(3,3);
M(5,5) = 4/3*pi*R3^3*rho_metal;
M(6,6) = M(5,5);

% Viscous damping matrix
C = zeros(6,6);
C1 = 6*pi*visc*R1;
C2 = 6*pi*visc*R2;
C3 = 6*pi*visc*R3;
C(1,1) = C1;
C(2,2) = C1;
C(3,3) = C2;
C(4,4) = C2;
C(5,5) = C3;
C(6,6) = C3;

% Gravity
W = zeros(ndof,1);
W(2) = -4/3*pi*R1^3*rho*g;
W(4) = -4/3*pi*R2^3*rho*g;
W(6) = -4/3*pi*R3^3*rho*g;

% Initial DoF vector
q0 = zeros(ndof,1); % initial position
for c=1:N
    q0 ( 2*c - 1 ) = nodes(c,1); % x coordinate
    q0 ( 2*c ) = nodes(c,2); % y coordinate
end

u0 = zeros(ndof,1); % initial velocity

% Tolerance
tol = EI / RodLength^2 * 1e-3;

% Method Choice: Implicit or Explicit; Comment the other command!
method = "Implicit";
% method = "Explicit";

% Time step size
if method == "Implicit"
    dt = 1e-2; % 1e-1 for Q4, default: 1e-2
else
    dt = 1e-5; % 1e-4 for Q4, default: 1e-5
end

% Number of time steps y-velocity of the middle node
Nsteps = round( totalTime / dt ) + 1;
all_mid_y = zeros( Nsteps, 1); % y-position of R2
all_mid_v = zeros( Nsteps, 1); % y-velocity of R2
all_mid_y(1) = q0(4);
all_mid_v(1) = u0(4);

% Store information
Q = zeros(ndof, Nsteps);
Q(:,1) = q0;
U = zeros(ndof, Nsteps);
U(:,1) = u0;

% Time marching scheme
for c=2:Nsteps
 
%     fprintf('Time = %f\n', (c-1) * dt );
    
    if(method == "Implicit")

        %% Implicit Method
        q = q0; % Guess
        % Newton Raphson
        err = 10 * tol;
        while err > tol
            %% Inertia
            f = M / dt * ( (q-q0) / dt - u0 );
            J = M / dt^2;

            %% Elastic forces
            % Linear spring 1 between nodes 1 and 2
            xk = q(1);
            yk = q(2);
            xkp1 = q(3);
            ykp1 = q(4);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            f(1:4) = f(1:4) + dF;
            J(1:4,1:4) = J(1:4,1:4) + dJ;
            
            % Linear spring 2 between nodes 2 and 3
            xk = q(3);
            yk = q(4);
            xkp1 = q(5);
            ykp1 = q(6);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
            f(3:6) = f(3:6) + dF;
            J(3:6,3:6) = J(3:6,3:6) + dJ;
            
            % Bending spring at node 2
            xkm1 = q(1);
            ykm1 = q(2);
            xk = q(3);
            yk = q(4);
            xkp1 = q(5);
            ykp1 = q(6);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            f(1:6) = f(1:6) + dF;
            J(1:6,1:6) = J(1:6,1:6) + dJ;
            
            %% Viscous force
            f = f + C * ( q - q0 ) / dt;
            J = J + C / dt;
            
            %% Weight
            f = f - W;
            
            % Update
            q = q - J \ f;            
            err = sum(abs(f));
        end
        % Update
        u = (q - q0) / dt; % Velocity

    else
        %% Explicit Method
        q = q0; % Last step
        E_q = zeros(ndof,1);

        % Linear spring 1 between nodes 1 and 2
        xk = q(1);
        yk = q(2);
        xkp1 = q(3);
        ykp1 = q(4);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        E_q(1:4) = E_q(1:4) + dF;

        % Linear spring 2 between nodes 2 and 3
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
        E_q(3:6) = E_q(3:6) + dF;
        
        % Bending spring at node 2
        xkm1 = q(1);
        ykm1 = q(2);
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        curvature0 = 0;
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
        E_q(1:6) = E_q(1:6) + dF;

        % Update
        u = (M \ (W - E_q - C*u0 ) * dt + u0);
        q = u * dt + q0;
    end

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
    all_mid_y(c) = q(4);
    all_mid_v(c) = u(4);
end

%% Plotting
timeArray = (0:Nsteps-1) * dt; % [0 10] seconds

% Q1:
figure(2); % Shapes of structure
grid on;
hold on;
time_sep = [0, 0.01, 0.05, 0.10, 1.0, 10.0];
for i = 1:length(time_sep)
    ind = find(timeArray == time_sep(i));
    plot( Q(1:2:end,ind), Q(2:2:end,ind), 'o-','LineWidth',1.5); % plot 3 nodes
end
axis equal;
legend('t = 0','t = 0.01','t = 0.05','t = 0.10','t = 1.0','t = 10.0')
xlabel('x (meter)');
ylabel('y (meter)');
% title('Shapes of structure evolving');

figure(3); % Position and Velocity of R2
subplot(1,2,1)
plot(timeArray, Q(4,:), 'r-','LineWidth',1.5);
xlabel('t (sec)');
ylabel('y (meter)');
% title('Position of R_2');
grid on
subplot(1,2,2)
plot(timeArray, U(4,:), 'b-','LineWidth',1.5);
xlabel('t (sec)');
ylabel('v [meter/sec]');
% title('Velocity (along y-axis) of R_2');
grid on

% Q2:
figure(4) % Terminal velocity of the system
plot(Q(1:2:end,end), U(2:2:end,end), 'r*-','LineWidth',1.5);
grid on
ylim([-9 0]*1e-3)
xlabel('x (meter)');
ylabel('v [meter/sec]');
% title('Terminal velocity along y-axis');

% Q3 & Q4:
figure(5) % See what happens to velocities when R1 = R2 = R3 or dt changes
plot(timeArray,U(2,:),timeArray,U(4,:),timeArray,U(6,:),'LineWidth',1.5);
grid on
legend('R1','R2','R3');
xlabel('t (sec)');
ylabel('v [meter/sec]');
% title('Velocities of the system');


 
