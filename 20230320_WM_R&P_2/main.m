% subjID: e.g. S001
clear all;
clc;


SubjID   = 'test';
Ori_reward = 2;
sess_num = 1;
run_num  = 1;    

%%  
offset = [0,0];
CurrDir = pwd;
SetupRand;
% set_test_gamma;
HideCursor; 
parameters;

%%
% sample;
practice;
% wm_test;
% wm_train;

delete *.asv
