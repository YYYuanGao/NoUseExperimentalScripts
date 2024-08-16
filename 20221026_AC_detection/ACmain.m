clear all;
clc;

addpath(genpath('Palamedes1_10_11'));
addpath(genpath('set_gamma'));

%% subInfo
SubjID   = 'ysp';
SessName = 'preTT';      %Prac, preTT, postT
Sess_Num = 1; 
Run_Num  = 1; 
location_used = 4; 
% ori_used = 1;
sample_ratio = 0.3;

%%
offset = [0,0];
CurrDir = pwd;
SetupRand; 
% set_test_gamma;
HideCursor;
parameters;

%% 
sample; 
% prac;
% train;

delete *.asv
