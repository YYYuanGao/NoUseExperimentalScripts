clear all;
clc;

addpath(genpath('Palamedes1_10_11'));
addpath(genpath('set_gamma'));

%% subInfo
SubjID   = 'test';
SessName = 'pretest';
Sess_Num = 1; 
Run_Num  = 2; 
location_used = 4; 

%%
sp = BioSemiSerialPort();
% sp.sendTrigger()
offset = [0,0];
CurrDir = pwd;
SetupRand; 
set_test_gamma;
HideCursor;
parameters;

%% 
% sample; 
train;

delete *.asv
