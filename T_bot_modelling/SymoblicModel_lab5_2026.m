clearvars
close all
clc

%% 1. Define symbolic parameters

% Geometry
syms R0 r rp real

% Masses
syms M m mp Mc Mb real

% Gravity
syms g real

% Inertias
syms J J1 J2 J3 J4 J5 real

% Stiffness coefficients
syms k1 k2 k3 k4 k5 k6 k7 k8 real

% Damping coefficients
syms Bc Bb B0 B1 B2 B3 B4 B5 real

% Motor/electrical parameters
syms Kt Ke R real

% External forces/torques
syms Fc Ft Fs N real
syms TM1 TM2 Via Vib U E real

%% 2. Generalized coordinates
% State mapping (for lab5_2026.m):
%   x1=x,  x2=dx,  x3=z,  x4=dz
%   x5=ha, x6=dha, x7=hb, x8=dhb
%   x9=h1, x10=dh1, x11=h2, x12=dh2
%   x13=h3,x14=dh3, x15=h4, x16=dh4
%   x17=h5,x18=dh5, x19=g

syms x z ha hb h1 h2 h3 h4 h5 real
q = [x; z; ha; hb; h1; h2; h3; h4; h5];

syms dx dz dha dhb dh1 dh2 dh3 dh4 dh5 real
dq = [dx; dz; dha; dhb; dh1; dh2; dh3; dh4; dh5];

syms ddx ddz ddha ddhb ddh1 ddh2 ddh3 ddh4 ddh5 real
ddq = [ddx; ddz; ddha; ddhb; ddh1; ddh2; ddh3; ddh4; ddh5];

%% 3. Virtual displacements
syms sx sz sha shb sh1 sh2 sh3 sh4 sh5 real
delta_q = [sx; sz; sha; shb; sh1; sh2; sh3; sh4; sh5];

%% 4. Kinetic energy
% Each pulley (P1-P4) moves with the cart in x and rotates.
% Pulley P5 and center beam move in both x and z.

T1 = 1/2*Mb*(dz^2);                          % Center beam in z
T2 = 1/2*(Mb+Mc)*(dx^2);                     % Cart + center beam in x
T3 = 1/2*J*(dha^2);                          % Motor pulley A
T4 = 1/2*J*(dhb^2);                          % Motor pulley B
T5 = 1/2*J1*(dh1^2) + 1/2*m*(dx^2);         % Pulley P1 (rides on cart)
T6 = 1/2*J2*(dh2^2) + 1/2*m*(dx^2);         % Pulley P2
T7 = 1/2*J3*(dh3^2) + 1/2*m*(dx^2);         % Pulley P3
T8 = 1/2*J4*(dh4^2) + 1/2*m*(dx^2);         % Pulley P4
T9 = 1/2*J5*(dh5^2) + 1/2*mp*(dx^2+dz^2);   % Pulley P5 (top of beam)

T = T1+T2+T3+T4+T5+T6+T7+T8+T9;

%% 5. Potential energy
% 8 belt segments modelled as linear springs (see Fig. 3 of lab manual).

% Vertical belts (drive z via P3,P4)
V3 = 1/2*k3*(h3*r + z)^2;
V4 = 1/2*k4*(-h4*r + z)^2;

% Horizontal belts (drive x via motor pulleys A,B and P1-P4)
V1 = 1/2*k1*(h1*r + ha*R0 + x)^2;
V2 = 1/2*k2*(-h3*r - ha*R0 + x)^2;
V5 = 1/2*k5*(h4*r + hb*R0 - x)^2;
V6 = 1/2*k6*(-h2*r - hb*R0 - x)^2;

% Diagonal belts (connect P1,P2 to P5, drive z)
V7 = 1/2*k7*(h2*r + h5*rp - z)^2;
V8 = 1/2*k8*(-h1*r - h5*rp - z)^2;

% Gravity (only center beam and P5 change height)
V9  = -Mb*g*z;
V10 = -mp*g*z;

V = V1+V2+V3+V4+V5+V6+V7+V8+V9+V10;

%% 6. Lagrangian
L = T - V;

%% 7. Non-conservative generalized forces
% Virtual work includes:
%   - Nonlinear friction on cart (x) and center beam (z)
%   - Motor torques TM1 (motor A) and TM2 (motor B)
%   - Viscous damping on all pulleys

Qnc = sym(zeros(length(q),1));

sW = -(dx*Bc + Fc*tanh(dx*N) + Fs*exp(-Ft*dx.^2).*tanh(dx*N))*sx ...
     -(dz*Bc + Fc*tanh(dz*N) + Fs*exp(-Ft*dz.^2).*tanh(dz*N))*sz + ...
     (TM1 - B0*dha)*sha + (TM2 - B0*dhb)*shb ...
     - B1*dh1*sh1 - B2*dh2*sh2 - B3*dh3*sh3 - B4*dh4*sh4 - B5*dh5*sh5;

for i = 1:9
    ss = delta_q(i);
    Qnc(i) = simplify(diff(sW, ss));
end

disp('Non-conservative generalized forces Qnc:')
disp(Qnc)

%% 8. Lagrange equations
% d/dt(dL/ddq_j) - dL/dq_j = Qnc_j

% Partial derivative of L with respect to q (spring + gravity restoring forces)
dLdq = jacobian(L, q).';

% Partial derivative of L with respect to dq (momentum terms)
dLddq = jacobian(L, dq).';

% Time derivative of dL/ddq:
% Since L is quadratic in velocities, d/dt(dL/ddq) = M_inertia * ddq
% where M_inertia is the (constant) inertia matrix.
% jacobian(dLddq, dq) extracts M_inertia; multiply by ddq to get d/dt.
dt_dLddq = jacobian(dLddq, dq) * ddq;

disp('dLdq:');    disp(simplify(dLdq))
disp('dLddq:');   disp(simplify(dLddq))
disp('dt_dLddq:'); disp(simplify(dt_dLddq))

%% 9. Solve for accelerations (ddq)
% Lagrange equation: dt_dLddq - dLdq - Qnc = 0
% Rearranged:        M_inertia * ddq = dLdq + Qnc

EOM = simplify(dt_dLddq - dLdq - Qnc);

% Extract inertia matrix (coefficient of ddq vector)
M_inertia = jacobian(EOM, ddq);

% RHS: everything that does not involve accelerations
RHS = simplify(-(EOM - M_inertia*ddq));

% Solve: ddq = M_inertia \ RHS
ddq_sol = simplify(M_inertia \ RHS);

disp('=== EQUATIONS OF MOTION (Accelerations) ===')
fprintf('\nddx  (-> A2 row):\n'); disp(ddq_sol(1))
fprintf('\nddz  (-> A4 row):\n'); disp(ddq_sol(2))
fprintf('\nddha (-> A6 row):\n'); disp(ddq_sol(3))
fprintf('\nddhb (-> A8 row):\n'); disp(ddq_sol(4))
fprintf('\nddh1 (-> A10 row):\n'); disp(ddq_sol(5))
fprintf('\nddh2 (-> A12 row):\n'); disp(ddq_sol(6))
fprintf('\nddh3 (-> A14 row):\n'); disp(ddq_sol(7))
fprintf('\nddh4 (-> A16 row):\n'); disp(ddq_sol(8))
fprintf('\nddh5 (-> A18 row):\n'); disp(ddq_sol(9))

%% 10. Linear state-space A-matrix extraction
% For the linear model: set Fc=Fs=0, substitute motor back-EMF
%   TM1 = -(Kt*Ke/R)*dha   (input Vin is captured in B matrix)
%   TM2 = -(Kt*Ke/R)*dhb

ddq_linear = subs(ddq_sol, [TM1, TM2, Fc, Fs], ...
                  [-(Kt*Ke/R)*dha, -(Kt*Ke/R)*dhb, 0, 0]);
ddq_linear = simplify(ddq_linear);

% Full state vector (19 states, g is treated as a constant state = 9.8)
state_vars = [x; dx; z; dz; ha; dha; hb; dhb; h1; dh1; ...
              h2; dh2; h3; dh3; h4; dh4; h5; dh5; g];

disp('=== A-MATRIX ROWS (linear, even-indexed rows for accelerations) ===')
row_labels = {'A2','A4','A6','A8','A10','A12','A14','A16','A18'};
for i = 1:9
    Arow = jacobian(ddq_linear(i), state_vars.');
    fprintf('\n%s =\n', row_labels{i});
    disp(simplify(Arow))
end