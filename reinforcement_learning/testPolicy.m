clear;
close all;

%% Input data:
% Fictional action:
actions = [-1,0,1];
% Fictional state:
n = 5;
states = nchoosek(1:n,n-1);
% Exploration rate:
epsilon = 0.8;

%% Test case: exact basis function
policy = Policy(actions,states,epsilon);