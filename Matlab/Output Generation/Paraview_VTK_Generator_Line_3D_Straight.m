%% Create vtk file from txt output of C++ code
clearvars; close all; clc;
% Warnings pop up with each iteration of RN_3D_Paraview due to vecnorm function
warning off 
addpath('Functions')
addpath('../Input File Generation/Functions')
addpath('../Formulation/Functions')
addpath('../Formulation/German Formulation/With Torsion/Functions_RN Parallel/3D_New')

%% Find Deformed shape
% Select deformation type for results
deformation_type = 'BIAX';
% Pick fiber number and Case number to analyze
fiber_num=313; Case=2; 

% Assign other parameters the same as in "Input_File_Generation.m"
power=4; ele_size=0.4/(2^(power));
order=5; end_cpt_factor=1/5; trim_factor=0.01;

% File/Folder names
MS=strcat('_MS',num2str(power)); MS_case=strcat('__MS',num2str(power));
folder_name=num2str(fiber_num);
network_name=strcat(folder_name,'_',deformation_type,'_Case',num2str(Case),'_P',num2str(order),MS_case);
disp(strcat('Visualizing: ',network_name))
path_c = '../../Results/Outputs/';

% Load network
load(['../Networks/Data/3D/Straight/',num2str(fiber_num),'/Case',num2str(Case)]);

% Load results file
MeshFileName=['Straight_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS];
ResultFileName=[MeshFileName,'_Result.txt'];
ResultFilePath_c=[path_c,ResultFileName];
CPTS=importdata(ResultFilePath_c);
filepath = ['../../Results/Paraview/3D/',num2str(fiber_num),'_',num2str(Case),'/',deformation_type];

if ~exist(filepath,'dir')
    mkdir(filepath);
end
vtk_file = [filepath,'/Straight_',num2str(fiber_num),'_',deformation_type];

%% Generate Mesh
[ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator_3D( Segment,Junction,ele_size,order,end_cpt_factor,trim_factor,'of' );
step=1;
%% Material Properties.
% Fiber radius
r=0.001;
% C/S area,Young's Modulus & Density
Area=pi*(r)^2; E=10^6; v=0.2;
G=0.5*(1+v)*E; I2=(pi/4)*r^4; I3=I2; Ip=I3+I2;
Mat = [Area E G I2 I3 Ip];

% P stores the coordinates of each control point in its reference state
P=CPTS(1:4,:);

[XP,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,P,Patch,A0,Mat,ele_size);
% Group each increment into an individual cell to avoid broadcast variable
% warning if using parfor
num_increments = size(CPTS,1)/4;
CPTS_cell = cell(num_increments);
for i=1:num_increments
    CPTS_cell{i}=CPTS(4*(i-1)+1:4*(i-1)+4,:);
end


% XP gives the coordinates of control points in their reference config
% XQ gives the coordaintes of control points in their current config
% In each loop, we update the value of XQ to get the coordinates in
% the ith increment, and use those to calculate the displacements
% off the fibers in Write_VTK_Network_Line.

% Option to use parfor - just uncomment the below lines to do so
% parpool('local',8);
% par
for i=1:num_increments
    Q = CPTS_cell{i};
    [XQ,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,Q,Patch,A0,Mat,ele_size);
    Write_VTK_Network_Line(XP,XQ,EPS,KAPPA,ME,BE,TE,SE,vtk_file,i);
end
% p = gcp;
% delete(gcp('nocreate'))