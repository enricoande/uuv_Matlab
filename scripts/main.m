% main.m     e.anderlini@ucl.ac.uk     07/11/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program runs reinforcement learning control for an UUV.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

%% Initialize UUV object:
load('rov.mat');
uuv = Uuv(rov);