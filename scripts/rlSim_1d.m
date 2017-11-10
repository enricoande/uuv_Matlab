% rlSim_1d.m     e.anderlini@ucl.ac.uk     09/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates the motions of an UUV under reinforcement learning
% control (in a single degree of freedom).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clear global;
close all;

%% Allocate global memory to speed up code:
global lspi samples;

%% Fix the seed to the random number generator:
rng(0);

%% Input data:
% Simulation data:
dt   = 0.01;              % time step (s)
tEnd = 20;                % simulation end time (s)
t    = (0:dt:tEnd);       % time stamp (s)
% Physical properties:
nu_c = [0;0;0;0;0;0];     % current velocity (m/s)
% Waypoints and circle of acceptance:
waypoints = [0,0,0];
r = 0.1;
% Desired depth value:
z_d = 1;                  % desired depth value (m)

%% Initialize UUV object:
load('rov.mat');
uuv = Uuv(rov,dt);

%% Initialize policy and LSPI objects:
% Input data:
epsilon = 0.8;            % exploration rate
exploration_decay = 0.99; % exploration rate decay
delta = 0.95;             % discount factor
n_episodes = 49;         % no. episodes

% Action space:
actions = ((-600:200:600)/60)';   % n (rps)
% State space:
centres = build_statelist((-5:0.5:5),(-2:0.2:2));
mu = 0.2;

% Initialize the policy object:
policy = Policy(actions,centres,epsilon,'basis_rbf',mu);
% Initialize the LSPI object:
lspi = Lspi(policy,delta,10,1e-5,1);
% Initialize the output variable:
duration = zeros(n_episodes,1);

%% Run LSPI for the UUV control:
tic;
for i=1:n_episodes
    % Run an episode:
    [x,f,dur] = episode_1d(uuv,t,0.1,nu_c,z_d);
    % Store episode duration:
    duration(i)= dur;
    
    % Update the policy:
    if mod(i,50) == 0
        % Keep only one copy of each row:
        samples = unique(samples,'rows');
        % Update the policy:
        lspi = lspi.lstdq(samples);
    end
    
    % Update the epxloration rate:
    lspi.policy = lspi.policy.update_epsilon(exploration_decay);
end
toc;

%% Motions post-processing:
% plotMotions(t,x);
% plotForces(t,f);
% plotPath(x,waypoints);
% animateAUV(t,x,50,1,8);

%% Reinforcement learning post-processing:
plotConvergence(duration);