%% This scripts introduces undulations based on chosen parameters in existing
% straight fiber networks and saves new networks.
clearvars; close all; clc;

addpath Functions
%% Code
% Specifiy network
num_fibers = 350;
Case=1;
New_Case=1;
% Row vector containing periods of sin undulations that can be introduced in the fiber.
period = [0.5,1,1.5,2];
num_periods=size(period,2);
% Row vector used to select periods in local y and z that a fiber can take
% based on its length. The longer the fiber, the more options it has.
select_period=[0,0.1,0.2,0.4]; 

network_path= ['../Networks/Data/3D/Undulated/',num2str(num_fibers),'/Case',num2str(Case)];
new_path= ['../Networks/Data/3D/Undulated/',num2str(num_fibers)];
if ~exist(new_path,'dir')
    mkdir(new_path)
end
new_network_path = [new_path,'/Case',num2str(Case),'_',num2str(New_Case)];
Undulate_3D_Network(network_path,new_network_path,period,select_period);

