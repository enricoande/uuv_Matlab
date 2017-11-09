% rovSim.m     e.anderlini@ucl.ac.uk     08/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program simulates the motions of an UUV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%% Input data:
% Simulation data:
dt   = 0.01;            % time step (s)
tEnd = 20;              % simulation end time (s)
t    = (0:dt:tEnd);     % time stamp (s)
nT   = length(t);       % no. time steps
% Physical properties:
nu_c = [0;0;0;0;0;0];   % current velocity (m/s)
n = [0;0;0;0;0];        % propellers' revolutions (rps)
% Waypoints and circle of acceptance:
waypoints = [0,0,0;
             2,0,0;
             2,4,0];
r = 0.1;

%% Initialize UUV object:
load('rov.mat');
uuv = Uuv(rov,dt);

%% Initialization:
x = zeros(12,nT);       % state vector
f = zeros(24,nT);       % force vector

%% Simulate the motions of the UUVs:
tic;
for i=2:nT
    [tmp_x,uuv] = uuv.update_dynamics(x(:,i-1),n,nu_c);
    x(:,i) = tmp_x;
    f(:,i) = [uuv.tau;uuv.f_h;uuv.f_d;uuv.f_c];
end
toc;
% Flip x and f to correct orientation:
x = x';
f = f';

%% Post-processing:
plotMotions(t,x);
plotForces(t,f);
plotPath(x,waypoints);
% animateAUV(t,x,50,1,8);