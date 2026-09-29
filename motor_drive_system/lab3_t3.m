%% SRT700 Lab 3 - Task 3: Motor Drive Two-Mass System
% State space model of motor driving two inertias via flexible shaft
% Parameters: J1=10/9, J2=10, c=0.1, k=1, kI=1

clear; close all; clc;

fprintf('=== TASK 3: Motor Drive Two-Mass System ===\n\n');

%% Parameters
J1 = 10/9;
J2 = 10;
c  = 0.1;
k  = 1;
kI = 1;

% Undamped natural frequency
w0 = sqrt(k*(J1+J2)/(J1*J2));
fprintf('Natural frequency w0 = sqrt(k*(J1+J2)/(J1*J2)) = %.6f rad/s\n\n', w0);

%% Build state-space matrix A (4x4)
% States: x1=phi1, x2=phi2, x3=w1/w0, x4=w2/w0
A = [0,              0,             w0,        0;
     0,              0,             0,         w0;
    -k/(w0*J1),  k/(w0*J1),   -c/J1,     c/J1;
     k/(w0*J2), -k/(w0*J2),    c/J2,    -c/J2];

% Input matrix B (control: motor current I)
B_I = [0; 0; kI/(w0*J1); 0];

% Disturbance matrix B_d (disturbance torque on second rotor)
B_d = [0; 0; 0; 1/(w0*J2)];

fprintf('A matrix:\n');
disp(A);
fprintf('B (control) vector:\n');
disp(B_I');
fprintf('B_disturbance vector:\n');
disp(B_d');

%% Part (a): Verify open-loop eigenvalues
ol_eigs = eig(A);
fprintf('Part (a): Open-loop eigenvalues:\n');
for i = 1:4
    fprintf('  lambda_%d = %.6f + %.6fi\n', i, real(ol_eigs(i)), imag(ol_eigs(i)));
end
fprintf('Expected: 0, 0, -0.05 +/- i\n\n');

%% Part (b): State feedback design
% Desired eigenvalues: -2, -1, -1+i, -1-i
desired_poles = [-2, -1, -1+1i, -1-1i];
fprintf('Part (b): Desired closed-loop poles:\n');
fprintf('  -2, -1, -1+i, -1-i\n\n');

% Check controllability
Co = ctrb(A, B_I);
rk = rank(Co);
fprintf('Controllability matrix rank: %d / 4  => System is %s\n\n', ...
    rk, ternary(rk==4,'CONTROLLABLE','NOT fully controllable'));

% Compute state feedback gain K using pole placement
K = place(A, B_I, desired_poles);
fprintf('State feedback gain K:\n');
fprintf('  K = [%.4f, %.4f, %.4f, %.4f]\n\n', K(1), K(2), K(3), K(4));

% Verify closed-loop eigenvalues
A_cl = A - B_I*K;
cl_eigs = eig(A_cl);
fprintf('Closed-loop eigenvalues (verification):\n');
for i = 1:4
    fprintf('  lambda_%d = %.4f + %.4fi\n', i, real(cl_eigs(i)), imag(cl_eigs(i)));
end
fprintf('\n');

%% Part (c): Simulation
% Output: y = phi2 = x2  (angle of second mass)
C_out = [0, 1, 0, 0];   % observe phi2
D_out = 0;

% Steady-state pre-compensator for reference tracking
% At DC: x_ss = -(A_cl)^-1 * B_I * N * r,  y_ss = C*x_ss = 1 => N = -1/(C*inv(A_cl)*B_I)
N = -1 / (C_out * (A_cl \ B_I));
fprintf('Pre-compensator gain N = %.4f\n\n', N);

t = linspace(0, 25, 5000);

%--- Simulation 1: Step in reference signal ---
% dx = A_cl*x + B_I*(N*r),  y = C*x
sys_ref = ss(A_cl, B_I*N, C_out, D_out);
[y_ref, t_ref, x_ref] = step(sys_ref, t);

%--- Simulation 2: Step disturbance on second rotor ---
% u = -Kx (no reference, just disturbance)
% dx = (A - B_I*K)*x + B_d*Td = A_cl*x + B_d*Td
sys_dist = ss(A_cl, B_d, C_out, D_out);
[y_dist, t_dist] = step(sys_dist, t);

%--- Plot ---
figure('Name','Task 3 - Motor Drive Closed Loop','Position',[100 100 1100 600]);

% Reference step response
subplot(2,2,[1 2]);
plot(t_ref, y_ref, 'b', 'LineWidth', 2); hold on;
plot(t_ref, x_ref(:,1), 'r--', 'LineWidth', 1.5);
yline(1, 'k:', 'Reference r=1', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('\phi_2 (rad)', 'FontSize', 12);
title('Task 3(c): Closed-Loop Response to Step Reference (r=1)', 'FontSize', 13);
legend('\phi_2 (output, 2nd mass angle)', '\phi_1 (1st mass angle)', ...
       'Location','best', 'FontSize', 10);
grid on;

% Disturbance response
subplot(2,2,[3 4]);
plot(t_dist, y_dist, 'm', 'LineWidth', 2);
yline(0, 'k:', 'Desired = 0', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 12);
ylabel('\phi_2 (rad)', 'FontSize', 12);
title('Task 3(c): Closed-Loop Response to Step Disturbance Torque T_d=1', 'FontSize', 13);
legend('Disturbance response \phi_2', 'Location','best', 'FontSize', 10);
grid on;

sgtitle('Task 3: Two-Mass Motor Drive with State Feedback', 'FontSize', 14);
saveas(gcf, 'Task3_motor_drive.png');

% Also plot all four states for reference step
figure('Name','Task 3 - All States','Position',[100 100 1000 600]);
state_labels = {'\phi_1 (rad)', '\phi_2 (rad)', '\omega_1/\omega_0', '\omega_2/\omega_0'};
for s = 1:4
    subplot(2,2,s);
    plot(t_ref, x_ref(:,s), 'LineWidth', 2);
    grid on;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel(state_labels{s}, 'FontSize', 10);
    title(['State x_' num2str(s) ': ' state_labels{s}], 'FontSize', 11);
end
sgtitle('Task 3: All States - Step Reference Response', 'FontSize', 13);
saveas(gcf, 'Task3_all_states.png');

fprintf('Task 3 complete.\n');

%% Helper
function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end