function [x,f,duration] = episode_1d(uuv,t,dt,nu_c,z_d)
% episode_1d.m     e.anderlini@ucl.ac.uk     09/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function runs an episode of LSPI for the control of the UUV in one
% degree of freedom.
% Input:
% uuv:   object of UUV dynamics;
% t:     time series;
% dt:    time step of reinforcement learning control;
% nu_c:  current velocity vector;
% z_d:   desired depth value. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Allocate global memory to speed up code:
global lspi samples;

%% Initialization:
nT   = length(t);                    % no. time steps
x = zeros(12,nT);                    % state vector
f = zeros(24,nT);                    % force vector
n = zeros(uuv.propulsion.n_prop,1);  % propellers' revolutions  
s = [];                              % state variable
start_time = [];                     % start time variable

%% Simulate the motions of the UUVs:
for i=2:nT
    % Store output variable:
    duration = t(i);
    
    % Update system dynamics and store force terms:
    [tmp_x,uuv] = uuv.update_dynamics(x(:,i-1),n,nu_c);
    x(:,i) = tmp_x;
    f(:,i) = [uuv.tau;uuv.f_h;uuv.f_d;uuv.f_c];
    
    % Apply action only every 0.1 s:
    if mod(t(i),dt)==0
        % Store the new state:
        sp = [x(3,i),x(9,i)];

        % Find reward:
        [r,flag] = get_reward_1d(sp,t(i),z_d);

        % Take a new action:
        ap = lspi.policy.select_action(sp);

        % Store sample:
        if ~isempty(s)
            sample = [s,a,r,sp];
            samples = [samples;sample];
        end

        % Update state & action:
        s = sp;
        a = ap;

        % Check if simulation can be stopped:
        if flag && isempty(start_time)
            start_time = t(i);
        end
        if ~isempty(start_time)
            if t(i)>start_count+2
                break;
            end
        end

        % Update propellers' revolutions:
        action = lspi.policy.actions(a);
        % Output desired values:
        n = [0;action;action;0;0];
    end
end

%% Prepare output data:
% Flip x and f to correct orientation:
x = x';
f = f';

end