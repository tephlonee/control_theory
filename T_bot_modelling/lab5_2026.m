clear all;
close all;
clc;

%% =========================================================
%% PARAMETERS
%% =========================================================
r   = 0.0194;    % Pulley P1-P4 radius [m]
rp  = 0.0194;    % Pulley P5 radius [m]
R0  = 0.0194;    % Motor pulley radius [m]
M   = 0.001;     % Motor pulley mass [kg]
m   = 0.001;     % Pulley P1-P4 mass [kg]
mp  = 0.001;     % Pulley P5 mass [kg]
Mb  = 0.446;     % Center beam mass [kg]
Mc  = 4.07;      % Cart mass [kg]
g   = 9.8;       % Gravitational acceleration [m/s^2]
J   = 5.6e-5;    % Motor pulley inertia [kg.m^2]
J1  = 5.6e-5;    % Pulley P1 inertia
J2  = 5.6e-5;    % Pulley P2 inertia
J3  = 5.6e-5;    % Pulley P3 inertia
J4  = 5.6e-5;    % Pulley P4 inertia
J5  = 5.6e-5;    % Pulley P5 inertia

Bc  = 23.2224;   % Cart/beam viscous friction coefficient [Ns/m]
% Bc = 40;       % (alternative value)
Bb  = 46.4449;   % (unused - legacy)
B1  = 0.0021;    % Pulley P1 damping [Ns.m]
B2  = 0.0021;
B3  = 0.0021;
B4  = 0.0021;
B5  = 0.0021;
B0  = 0.0021;    % Motor pulley damping [Ns.m]

Fc  = 0.16;      % Coulomb friction coefficient (nonlinear model)
Fs  = 0.13;      % Stribeck friction coefficient
Ft  = 100000;    % Stribeck decay rate
N   = 1000;      % Smoothing factor for tanh approximation

Kt  = 0.023;     % Motor torque constant [Nm/A]
Ke  = 0.023;     % Back-EMF constant [V/(rad/s)]
R   = 0.71;      % Armature resistance [Ohm]

k1  = 17810;     % Belt segment stiffnesses [N/m]
k2  = 17810;
k5  = 17810;
k6  = 17810;
k3  = 12565;
k4  = 12565;
k7  = 25013;
k8  = 25013;

%% =========================================================
%% DERIVED QUANTITIES
%% =========================================================

% Total effective masses
Mx = Mb + Mc + 4*m + mp;   % Mass in x-direction
Mz = Mb + mp;               % Mass in z-direction

% Back-EMF + bearing damping for motor pulleys
Beff_motor = B0 + Kt*Ke/R;  % Total effective damping on motor shaft

%% =========================================================
%% STATE VECTOR (19 states)
%% =========================================================
% x  = [x, dx, z, dz, ha, dha, hb, dhb, h1, dh1,
%        h2, dh2, h3, dh3, h4, dh4, h5, dh5, g]
% idx=  1   2   3   4   5   6   7   8   9  10
%       11  12  13  14  15  16  17  18  19
%
% Notation: ha=phi_A, hb=phi_B, h1..h5=phi_1..phi_5

%% =========================================================
%% A MATRIX ASSEMBLY
%% =========================================================
% Odd rows (1,3,5,...,17): dx_j/dt = velocity state  (trivial)
% Even rows (2,4,6,...,18): d(vel)/dt = acceleration (from Lagrange EOM)
% Row 19: dg/dt = 0 (g is a constant parameter treated as a state)
%
% Equations of motion derived from Lagrangian (see SymbolicModel_lab5_2026.m):
%
% ddx  = [-(k1+k2+k5+k6)*x + R0*(k2-k1)*ha + R0*(k5-k6)*hb
%          - r*(k1*h1 + k6*h2 - k2*h3 - k5*h4) - Bc*dx] / Mx
%
% ddz  = [-(k3+k4+k7+k8)*z - r*(k8*h1 - k7*h2 + k3*h3 - k4*h4)
%          + rp*(k7-k8)*h5 + (Mb+mp)*g - Bc*dz] / Mz
%
% ddha = [-(k1+k2)*R0^2*ha - R0*r*(k1*h1 + k2*h3) + R0*(k2-k1)*x
%          - Beff_motor*dha] / J   [+Kt/R * Va enters via B matrix]
%
% ddhb = [-(k5+k6)*R0^2*hb - R0*r*(k6*h2 + k5*h4) + R0*(k5-k6)*x
%          - Beff_motor*dhb] / J   [+Kt/R * Vb enters via B matrix]
%
% ddh1 = [-k1*R0*r*ha - (k1+k8)*r^2*h1 - k8*r*rp*h5 - k1*r*x - k8*r*z
%          - B1*dh1] / J1
%
% ddh2 = [-k6*R0*r*hb - (k6+k7)*r^2*h2 - k7*r*rp*h5 - k6*r*x + k7*r*z
%          - B2*dh2] / J2
%
% ddh3 = [-k2*R0*r*ha - (k2+k3)*r^2*h3 + k2*r*x - k3*r*z
%          - B3*dh3] / J3
%
% ddh4 = [-k5*R0*r*hb - (k4+k5)*r^2*h4 + k5*r*x + k4*r*z
%          - B4*dh4] / J4
%
% ddh5 = [-(k8*r*rp*h1 + k7*r*rp*h2) - (k7+k8)*rp^2*h5 + (k7-k8)*rp*z
%          - B5*dh5] / J5

% --- Row A1: dx1/dt = x2 ---
A1 = zeros(1,19);
A1(2) = 1;

% --- Row A2: ddx ---
% ddx = f(x[1]=x, x[2]=dx, x[5]=ha, x[7]=hb, x[9]=h1, x[11]=h2, x[13]=h3, x[15]=h4)
A2 = zeros(1,19);
A2(1)  = -(k1+k2+k5+k6) / Mx;    % x coefficient
A2(2)  = -Bc / Mx;                % dx coefficient (viscous friction)
A2(5)  = R0*(k2-k1) / Mx;        % ha coefficient (with k1=k2 this is 0)
A2(7)  = R0*(k5-k6) / Mx;        % hb coefficient (with k5=k6 this is 0)
A2(9)  = -k1*r / Mx;             % h1 coefficient
A2(11) = -k6*r / Mx;             % h2 coefficient
A2(13) =  k2*r / Mx;             % h3 coefficient
A2(15) =  k5*r / Mx;             % h4 coefficient
% A2(19) = 0  (no gravity in x)

% --- Row A3: dx3/dt = x4 ---
A3 = zeros(1,19);
A3(4) = 1;

% --- Row A4: ddz ---
% ddz = f(x[3]=z, x[4]=dz, x[9]=h1, x[11]=h2, x[13]=h3, x[15]=h4, x[17]=h5, x[19]=g)
A4 = zeros(1,19);
A4(3)  = -(k3+k4+k7+k8) / Mz;   % z coefficient
A4(4)  = -Bc / Mz;               % dz coefficient (viscous friction)
A4(9)  = -k8*r / Mz;             % h1 coefficient
A4(11) =  k7*r / Mz;             % h2 coefficient
A4(13) = -k3*r / Mz;             % h3 coefficient
A4(15) =  k4*r / Mz;             % h4 coefficient
A4(17) = (k7-k8)*rp / Mz;       % h5 coefficient (=0 since k7=k8)
A4(19) = 1;                      % g coefficient: (Mb+mp)*g / (Mb+mp) = g

% --- Row A5: dx5/dt = x6 ---
A5 = zeros(1,19);
A5(6) = 1;

% --- Row A6: ddha ---
% ddha = f(x[1]=x, x[5]=ha, x[6]=dha, x[9]=h1, x[13]=h3)
% Note: Beff_motor = B0 + Kt*Ke/R (back-EMF absorbed into A)
A6 = zeros(1,19);
A6(1)  = R0*(k2-k1) / J;         % x coefficient (=0 since k1=k2)
A6(5)  = -(k1+k2)*R0^2 / J;     % ha coefficient
A6(6)  = -Beff_motor / J;        % dha coefficient
A6(9)  = -R0*k1*r / J;           % h1 coefficient
A6(13) = -R0*k2*r / J;           % h3 coefficient

% --- Row A7: dx7/dt = x8 ---
A7 = zeros(1,19);
A7(8) = 1;

% --- Row A8: ddhb ---
% ddhb = f(x[1]=x, x[7]=hb, x[8]=dhb, x[11]=h2, x[15]=h4)
A8 = zeros(1,19);
A8(1)  = R0*(k5-k6) / J;         % x coefficient (=0 since k5=k6)
A8(7)  = -(k5+k6)*R0^2 / J;     % hb coefficient
A8(8)  = -Beff_motor / J;        % dhb coefficient
A8(11) = -R0*k6*r / J;           % h2 coefficient
A8(15) = -R0*k5*r / J;           % h4 coefficient

% --- Row A9: dx9/dt = x10 ---
A9 = zeros(1,19);
A9(10) = 1;

% --- Row A10: ddh1 ---
% ddh1 = f(x[1]=x, x[3]=z, x[5]=ha, x[9]=h1, x[10]=dh1, x[17]=h5)
A10 = zeros(1,19);
A10(1)  = -k1*r / J1;            % x coefficient
A10(3)  = -k8*r / J1;            % z coefficient
A10(5)  = -R0*k1*r / J1;        % ha coefficient
A10(9)  = -(k1+k8)*r^2 / J1;   % h1 coefficient
A10(10) = -B1 / J1;              % dh1 coefficient
A10(17) = -k8*r*rp / J1;        % h5 coefficient

% --- Row A11: dx11/dt = x12 ---
A11 = zeros(1,19);
A11(12) = 1;

% --- Row A12: ddh2 ---
% ddh2 = f(x[1]=x, x[3]=z, x[7]=hb, x[11]=h2, x[12]=dh2, x[17]=h5)
A12 = zeros(1,19);
A12(1)  = -k6*r / J2;            % x coefficient
A12(3)  =  k7*r / J2;            % z coefficient
A12(7)  = -R0*k6*r / J2;        % hb coefficient
A12(11) = -(k6+k7)*r^2 / J2;   % h2 coefficient
A12(12) = -B2 / J2;              % dh2 coefficient
A12(17) = -k7*r*rp / J2;        % h5 coefficient

% --- Row A13: dx13/dt = x14 ---
A13 = zeros(1,19);
A13(14) = 1;

% --- Row A14: ddh3 ---
% ddh3 = f(x[1]=x, x[3]=z, x[5]=ha, x[13]=h3, x[14]=dh3)
A14 = zeros(1,19);
A14(1)  =  k2*r / J3;            % x coefficient
A14(3)  = -k3*r / J3;            % z coefficient
A14(5)  = -R0*k2*r / J3;        % ha coefficient
A14(13) = -(k2+k3)*r^2 / J3;   % h3 coefficient
A14(14) = -B3 / J3;              % dh3 coefficient

% --- Row A15: dx15/dt = x16 ---
A15 = zeros(1,19);
A15(16) = 1;

% --- Row A16: ddh4 ---
% ddh4 = f(x[1]=x, x[3]=z, x[7]=hb, x[15]=h4, x[16]=dh4)
A16 = zeros(1,19);
A16(1)  =  k5*r / J4;            % x coefficient
A16(3)  =  k4*r / J4;            % z coefficient
A16(7)  = -R0*k5*r / J4;        % hb coefficient
A16(15) = -(k4+k5)*r^2 / J4;   % h4 coefficient
A16(16) = -B4 / J4;              % dh4 coefficient

% --- Row A17: dx17/dt = x18 ---
A17 = zeros(1,19);
A17(18) = 1;

% --- Row A18: ddh5 ---
% ddh5 = f(x[3]=z, x[9]=h1, x[11]=h2, x[17]=h5, x[18]=dh5)
% Note: no x dependence, no motor (P5 is passive)
A18 = zeros(1,19);
A18(3)  = (k7-k8)*rp / J5;      % z coefficient (=0 since k7=k8)
A18(9)  = -k8*r*rp / J5;        % h1 coefficient
A18(11) = -k7*r*rp / J5;        % h2 coefficient
A18(17) = -(k7+k8)*rp^2 / J5;  % h5 coefficient
A18(18) = -B5 / J5;              % dh5 coefficient

% --- Row A19: dg/dt = 0 ---
A19 = zeros(1,19);

%% Assemble full A matrix
A = [A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;A13;A14;A15;A16;A17;A18;A19];

%% =========================================================
%% B MATRIX (input voltage to motor accelerations)
%% =========================================================
% Inputs: u = [Va; Vb]  (voltage to motor A and motor B)
% Motor torque: TM = (Kt/R)*V  (at rest, no back-EMF; back-EMF is in A)
% Effect: ddha = (Kt/(J*R))*Va  => B(6,1)
%         ddhb = (Kt/(J*R))*Vb  => B(8,2)
B = zeros(19,2);
B(6,1) = Kt/(J*R);    % Motor A voltage -> angular acceleration of ha
B(8,2) = Kt/(J*R);    % Motor B voltage -> angular acceleration of hb

%% =========================================================
%% C MATRIX (output selection)
%% =========================================================
% Two output configurations:
%   CC : output x (state 1) and z (state 3) - Cartesian positions
%   CC1: output ha (state 5) and hb (state 7) - motor encoder angles (physically measured)

CC = zeros(2,19);
CC(1,1) = 1;    % Output 1: x position
CC(2,3) = 1;    % Output 2: z position

CC1 = zeros(2,19);
CC1(1,5) = 1;   % Output 1: motor A angle (encoder reading)
CC1(2,7) = 1;   % Output 2: motor B angle (encoder reading)

D = zeros(2,2);

%% =========================================================
%% RELATIVE GAIN ARRAY (RGA) - coupling analysis
%% =========================================================
% RGA at DC (omega=0) tells us how coupled the x and z channels are.
% RGA(i,j) = G(i,j)*[G^-1(j,i)]: diagonal dominance -> easier to control.
% For the 18-state linear system (excluding g row/column):

g0 = ss(A(1:18,1:18), B(1:18,:), CC(:,1:18), []);
G  = freqresp(g0, 0);        % DC gain matrix (2x2)
RGA = G .* pinv(G).'         % Relative Gain Array

disp('RGA at DC:')
disp(RGA)
disp('RGA diagonal values close to 1 indicate low coupling.')

%% =========================================================
%% GRAVITY COMPENSATION VOLTAGE U0
%% =========================================================
% The center beam (mass Mb) and pulley P5 (mass mp) hang in the z-direction.
% Their combined weight must be balanced by motor torques through the belt.
%
% Force balance at static equilibrium (z=0 desired):
%   - Total weight: W = (Mb + mp)*g
%   - Each motor supports W/2 through the belt
%   - Belt force -> motor torque: TM = (W/2) * R0
%     (force at pulley radius R0 provides the required tension)
%   - Motor voltage: TM = (Kt/R)*U0  =>  U0 = TM*R/Kt
%
% This U0 is added as a constant feedforward to BOTH motors.

W   = (Mb + mp) * g;                 % Total weight to support [N]
TM0 = (W / 2) * R0;                  % Required torque per motor [Nm]
U0  = TM0 * R / Kt;                  % Required voltage per motor [V]

fprintf('\n--- Gravity Compensation ---\n')
fprintf('Weight of beam + P5:  W  = %.4f N\n', W)
fprintf('Torque per motor:     TM0 = %.6f Nm\n', TM0)
fprintf('Compensation voltage: U0  = %.4f V\n', U0)
fprintf('Apply U0 = %.4f V to BOTH motors to hold beam at z = 0\n\n', U0)

%% =========================================================
%% PRINT A MATRIX SUMMARY
%% =========================================================
fprintf('A matrix assembled (%dx%d).\n', size(A,1), size(A,2))
fprintf('Non-zero entries in A: %d\n', nnz(A))
fprintf('Eigenvalues of A(1:18,1:18) (open-loop poles):\n')
ev = eig(A(1:18,1:18));
disp(ev)

%% =========================================================
%% CREATE STATE-SPACE OBJECTS FOR SIMULATION
%% =========================================================
% Linear state-space (full 18-state, g treated as constant input via x19=9.8)
% For simulation, set initial condition x0(19)=g=9.8

% sys_lin: linear system (x and z outputs)
sys_lin = ss(A(1:18,1:18), [B(1:18,:), A(1:18,19)], CC(:,1:18), zeros(2,3));
% The third input column represents the gravity input (constant 9.8)

% sys_lin_encoder: same but encoder (motor angle) outputs
sys_lin_enc = ss(A(1:18,1:18), [B(1:18,:), A(1:18,19)], CC1(:,1:18), zeros(2,3));

fprintf('\nLinear state-space model created.\n')
fprintf('Inputs: [Va, Vb, g_input]\n')
fprintf('Outputs (sys_lin): [x, z]\n')
fprintf('Outputs (sys_lin_enc): [ha, hb]\n')