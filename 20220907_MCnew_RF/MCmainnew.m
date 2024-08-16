clear all;
clc;

addpath(genpath('Palamedes1_10_11'));
addpath(genpath('set_gamma'));

%% subInfo
SubjID   = '111';
SessName = 'Train';
Sess_Num = 1; 
Run_Num  = 1; 
location_used = 1; 
ori_used = 2;

% %% priors
% prior_range =  0:0.5:8;  %-3:.01:3
% prior_sd = PAL_pdfNormal(Param.RF.alphas,4,2);

%%
offset = [0,0];
CurrDir = pwd;
SetupRand; 
% set_test_gamma;
HideCursor;
parameters;

%% 
sample; 
% MCtrain;

delete *.asv
