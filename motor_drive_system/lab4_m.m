
%% System Parameter Definition
J1 = 10 / 9;
J2 = 10;
c  = 0.1;
k  = 1;
kI = 1;

% Calculate the scaling frequency omega_0
omega0 = k * (J1 + J2) / (J1 * J2);

%% Task 1: Building the Continuous-Time State-Space Matrices

% Matrix A (4x4)
% Row 1: dx1/dt = omega0 * x3
% Row 2: dx2/dt = omega0 * x4
% Row 3: dx3/dt equations (normalized by omega0)
% Row 4: dx4/dt equations (normalized by omega0)
A = [          0,            0,      omega0,           0;
    0,            0,           0,      omega0;
    -k/(omega0*J1),  k/(omega0*J1),  -c/J1,        c/J1;
    k/(omega0*J2), -k/(omega0*J2),   c/J2,       -c/J2];

% Matrix B (4x1) - Input current (I) path
B = [0; 
    0; 
    kI / (omega0 * J1); 
    0];

% Matrix Bd (4x1) - Load disturbance (Td) path
% Crucial for Task 4, 6, and 7 simulation!
Bd = [0; 
    0; 
    0; 
    1 / (omega0 * J2)];

% Combined input matrix for standard Simulink block parameters: [B, Bd]
B_plant = [B, Bd];

% Matrix C (1x4) - Output is phi_2 (x2)
C = [0, 1, 0, 0];

% Matrix D (1x2) - Direct feedthrough for [Input, Disturbance]
D = [0, 0];

%% Display Verification
disp('--- Matrix A ---'); disp(A);
disp('--- Matrix B (Control Input) ---'); disp(B);
disp('--- Matrix Bd (Disturbance Input) ---'); disp(Bd);

%% %% TASK 1: State Feedback Pole Placement
poles_1_1 = [-0.8, -0.7, -0.6, -0.5];
K_1_1     = place(A, B, poles_1_1);

poles_1_2 = [-1.3, -1.2, -1.1, -1.0];
K_1_2     = place(A, B, poles_1_2);

% We use the fast poles configuration for the subsequent steps
K = K_1_2; 

%% %% TASK 2 & 3: Reference Feedforward Gain
% Calculate kr to eliminate steady-state error when Td = 0
kr = 1 / (-C * inv(A - B*K) * B);

%% %% TASK 4: Augmented Integral Control System
% Augmented matrices incorporating state x5 where dot(x5) = r - y
A_aug = [A,         zeros(4,1);
    -C,        0];
B_aug = [B; 0];

poles_aug = [-1.3, -1.2, -1.1, -1.0, -1.5]; % Added an extra stable pole for x5
K_aug_all = place(A_aug, B_aug, poles_aug);

% Split into state feedback gain (K_i_reg) and integral gain (Ki)
K_i_reg = K_aug_all(1:4); 
Ki      = -K_aug_all(5); % Note sign depending on block diagram mapping

%% %% TASK 5 & 6: State Observer Gain
poles_obs_slow = [-2.3, -2.2, -2.1, -2.0];
L_slow         = place(A', C', poles_obs_slow)';

poles_obs_fast = [-3.3, -3.2, -3.1, -3.0];
L_fast         = place(A', C', poles_obs_fast)';

% Default observer gain configuration for final loops
L = L_fast; 

%% %% Display Variables in Command Window
fprintf('--- Lab Parameters Successfully Calculated ---\n');
disp('State Feedback Gain K (Task 1.2/2/3):'), disp(K);
disp('Feedforward Gain kr (Task 2/3):'), disp(kr);
disp('Integral Controller Feedback Gains (Task 4):');
fprintf('  State component (K_i_reg): %s\n', mat2str(K_i_reg, 4));
fprintf('  Integral component (Ki):    %.4f\n', Ki);
disp('Observer Gain L (Task 5.2/6/7):'), disp(L);