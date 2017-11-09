function [r,stop] = get_reward_1d(s,duration,desired_depth)
% get_reward_1d.m      e.anderlini@ucl.ac.uk     05/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns the reward for the reinforcement learning problem
% of the control of the UUV in one dimension.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize reward & stop condition:
r = 0;
stop = false;

if abs(s(1)-desired_depth)<0.01
    r = r + 100;
    if abs(s(2))<0.01
        r = r + 1e03/duration;
        stop = true;
    end
end

end