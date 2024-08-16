 clear;
clc; 

addpath(genpath('Palamedes1_10_11'));                           x
addpath(genpath('set_gamma'));

%% subInfo
SubjID   = 'gy';
Sess_Num = 1;
Run_Num  = 1;  
%%
offset = [0,0];
CurrDir = pwd;
SetupRand;
% set_test_gamma;
HideCursor;
parameters;

%% 
% sample;
% sample2;
train; 
%  test;
delete *.asv
