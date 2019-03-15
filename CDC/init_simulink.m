%% Model parameters of the crazyflie 2.0
% Parameters
g = 9.81;       % m/s^2
m = 0.0331;      % kg
l = 0.046;      % m
k = 2.2e-8;  % N m s^2
b = 2e-9;       % N s^2
Im = 3e-6;    % kg*m^2
% I = [1.66e-5;    % kg*m^2
%      1.66e-5;    % kg*m^2
%      2.93e-5];   % kg*m^2
I = [2.3951e-5;    % kg*m^2
     2.3951e-5;    % kg*m^2
     3.2347e-5];   % kg*m^2
D = [0.92e-6;      % kg/s
	 0.91e-6;      % kg/s
	 1.03e-6];     % kg/s
 
Ct = 3.1582e-10;
Cd = 7.9379e-12;

% Thrust per rotor required to hover steadily
we = sqrt(m*g/(4*Ct));

% kalmanvar = [0.0005, 0.0005, 0.001234, 0.0014, 0.0014, 0.0015, 2e-5, 2e-5, 0.000256, 0.1, 0.1, 0.1];
% kalmanvar = [0.00382, 0.00382, 0.00749, 0.003, 0.003, 0.003, 2e-5, 2e-5, 8.144e-5, 0.01, 0.01, 0.01];
kalmanvar = [3.82e-5, 3.82e-5, 7.49e-5, 8e-4, 8e-4, 8e-4, 2e-5, 2e-5, 8.144e-5, 2e-6, 2e-6, 2e-6];

% new data: 0.0012, 0.00382, 0.00749, 0.00254, 0.00805, 0.00262, 2.477e-5,
% 1.265e-5, 8.144e-5

% variance measured on 3/16/2018:
% 2.0690e-04, 3.6354e-04, 1.9046e-3, 7.7743e-04, 7.8029e-04, 7.6584e-04, 
% 3.0958e-06, 4.1114e-06, 6.8328e-06, 1.1838e-06, 1.2549e-06, 1.1300e-06
%% Inner loop parameters for the crazyflie
% Loop rate at 500 Hz
inner_h = 0.002; 
controller_h = 0.004;
vicon_h = 0.01;
setpoint_h = 0.002;

% Omega saturation
inner_maxlim = 2500;
inner_minlim = 0;%500;

max_ang_ref = pi/8;

%% TDOA stuff
dwt_time = 499.2e6*128;

% anchorArray = [0.155, 0.190, 2.198;
%                4.500, 4.332, 0.195;
%                0.159, 0.780, 0.185;
%                4.498, 4.342, 2.185;
%                0.159, 4.365, 0.190;
%                4.498, 0.670, 0.190;
%                0.155, 4.245, 2.192;
%                4.495, 0.600, 2.195];

% IRL anchor positions
anchorArrayIRL = [8.25, 0.155, 3.312;
               0.155, 2.023, 3.208;
               9.827, 10.426, 3.265;
               0.155, 10.51, 3.2;
               10.086, 2.38, 0.212;
               0.155, 2.3, 0.211;
               9.827, 10.05, 0.212;
               0.155, 10.04, 0.214];
% offset = [4.0546, 5.4746, -0.2890];
anchorArray = anchorArrayIRL;
anchorArray(:,1) = anchorArray(:,1) - mean(anchorArray(:,1));
anchorArray(:,2) = anchorArray(:,2) - mean(anchorArray(:,2));
% anchorArray(:,3) = anchorArray(:,3) - mean(anchorArray(:,3));

% Initial states for the system
initial_xi = [0; 0; 0];
initial_xi2 = [-1; 0; 0];
initial_xi3 = [-2; 0; 0];
initial_xi4 = [-3; 0; 0];
initial_xi5 = [-4; 0; 0];
initial_xidot = [0; 0; 0];
initial_eta = [0; 0; 0];
initial_etadot = [0; 0; 0];

%% car stuff
d = 0.3;

max_acc = 100; %Have to find this
max_vel = 100;

max_steer_rate = 100;
max_steer_angle = pi/3;
init_theta = pi/6;