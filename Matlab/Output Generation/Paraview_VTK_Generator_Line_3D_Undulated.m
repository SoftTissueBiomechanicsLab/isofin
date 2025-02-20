%% Create vtk file from txt output of C++ code
clearvars; close all; clc;
% Warnings pop up with each iteration of RN_3D_Paraview due to vecnorm function
warning off 
addpath('Functions')
addpath('../Input File Generation/Functions')
addpath('../Formulation/Functions')
addpath('../Formulation/German Formulation/With Torsion/Functions_RN Parallel/3D_New')

%% Find Deformed shape
% Select deformation type and other parameters to select the desired
% results data. Be sure to use the same parameter values from your input
% file creation.
deformation_type = 'BIAX';
power=4; ele_size=0.4/(2^(power));
MS=strcat('_MS',num2str(power)); MS_case=strcat('__MS',num2str(power));
order=5; end_cpt_factor=1/5; trim_factor=0.01; Amp=0.05;
Case=1; New_Case=1; fiber_num=200;
folder_name=num2str(fiber_num);
network_name=strcat(folder_name,'_',deformation_type,'_Case',num2str(Case),'_P',num2str(order),MS_case);
disp(strcat('Visualizing: ',network_name))
path_c = '../../Results/Outputs/3D/';
%load network
load(['../Networks/Data/3D/Undulated/',num2str(fiber_num),'/Case',num2str(Case),'_',num2str(New_Case)]);

% % Desktop
MeshFileName=['Undulated_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_',num2str(New_Case),'_P',num2str(order),MS];
ResultFileName=[MeshFileName,'_Result.txt'];
ResultFilePath_c=[path_c,ResultFileName];
CPTS=importdata(ResultFilePath_c);

filepath = ['../../Results/Paraview/3D/',num2str(fiber_num),'/',deformation_type];
if ~exist(filepath,'dir')
    mkdir(filepath);
end
vtk_file = [filepath,'/Undulated_',num2str(fiber_num),'_',deformation_type];

%% Generate Mesh
[ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Undulated_Network_Mesh_Generator_3D( Segment,Segment_Period,Junction,Rotation,ele_size,order,period,Amp,trim_factor,'of' );


%% Material Properties.
step=1;
r=0.001;
Area=pi*(r)^2;E=10^6;v=0.2;G=0.5*(1+v)*E;I2=(pi/4)*r^4;I3=I2;Ip=I3+I2;% C/S area,Young's Modulus & Density
Mat(1)=Area;Mat(2)=E;Mat(3)=G;Mat(4)=I2;Mat(5)=I3;Mat(6)=Ip;
PtDensity=500;
P=CPTS(1:4,:);
[XP,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,P,Patch,A0,Mat,ele_size);

AtATime=20;
for j=1:ceil(((size(CPTS,1)/4))/AtATime)
    num=[AtATime*(j-1)+1,AtATime*j];
    parpool('local',6);
    for i=num(1):1:num(2)%(size(CPTS,1)/4)/2
        key=(i-1)*2+1;
        Q=CPTS(4*(key-1)+1:4*(key-1)+4,:);
        [XQ,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,Q,Patch,A0,Mat,ele_size);

        Write_VTK_Network_Line(XP,XQ,EPS,KAPPA,ME,BE,TE,SE,vtk_file,i);

    end
    p = gcp;delete(gcp('nocreate'));
end

%% Total energy
ME_network=zeros(1,size(CPTS,1)/4);BE_network=ME_network;TE_network=ME_network;SE_network=ME_network;
parfor i=1:size(CPTS,1)/4
    i
    Q=CPTS(4*(i-1)+1:4*(i-1)+4,:);
    [XQ,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,Q,Patch,A0,Mat,ele_size);
    t1=0;t2=0;t3=0;t4=0;
    for j=1:size(ME,2)
        t1=t1+sum(ME{j});t2=t2+sum(BE{j});
        t3=t3+sum(TE{j});t4=t4+sum(SE{j});
    end
    ME_network(i,1)=t1;BE_network(i,1)=t2;TE_network(i,1)=t3;SE_network(i,1)=t4;
end
n=size(BE_network,1);
strain=0.5*(0:1:n-1)./(n-1);
figure; hold on;
plot(strain,ME_network);plot(strain,BE_network);plot(strain,TE_network);plot(strain,SE_network);hold off;
figure; hold on;
plot(strain,ME_network./SE_network);plot(strain,BE_network./SE_network);plot(strain,TE_network./SE_network);hold off;