%% Calculate energy splits from txt output of C++ code - run after Plot_Results.m
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
% Select fiber number and case numbers for files
fiber_num=400; 
% 1 for straight, 2 for undulated network
bc_rotation = 1;
type = 1;
for Case= 1
    New_Case=1; 


% Assign other parameters the same as in "Input_File_Generation.m"
power=4; ele_size=0.4/(2^(power));
order=5; end_cpt_factor=1/5; trim_factor=0.01; Amp=0.025;

% File/Folder names
MS=strcat('_MS',num2str(power)); MS_case=strcat('__MS',num2str(power));
folder_name=num2str(fiber_num);
network_name=strcat(folder_name,'_',deformation_type,'_Case',num2str(Case),'_P',num2str(order),MS_case);
% path_c = '../../Results/Outputs/3D/';

switch type

    case 1    
    path_c = '../../Results/Outputs/'; %#ok<*UNRCH>
    load(['../Networks/Data/3D/Straight/',num2str(fiber_num),'/Case',num2str(Case)]);  
    % Load results file
    MeshFileName=['Straight_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS,'_',num2str(bc_rotation)];
    ResultFileName=[MeshFileName,'_Result.txt'];
    ResultFilePath_c=[path_c,ResultFileName];
    CPTS=importdata(ResultFilePath_c);
    
    filepath = ['../../Results/Paraview/3D/',num2str(fiber_num),'_',num2str(Case),'/',deformation_type];
    
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    vtk_file = [filepath,'/Straight_',num2str(fiber_num),'_',deformation_type];
    
    % Generate Mesh
    [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator_3D( Segment,Junction,ele_size,order,end_cpt_factor,trim_factor,'of' );

    case 2

    path_c = '../../Results/Outputs/'; 
    % Load network       
    load(['../Networks/Data/3D/Undulated/',num2str(fiber_num),'/Case',num2str(Case),'_',num2str(New_Case)]);
    
    
    % File paths and names
    MeshFileName=['Undulated_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_',num2str(New_Case),'_P',num2str(order),MS,'_',num2str(bc_rotation)];
    ResultFileName=[MeshFileName,'_Result.txt'];
    ResultFilePath_c=[path_c,ResultFileName];
    CPTS=importdata(ResultFilePath_c);
    filepath = ['../../Results/Paraview/3D/',num2str(fiber_num),'_',num2str(Case),'/',deformation_type];
    
    if ~exist(filepath,'dir')
        mkdir(filepath);
    end
    vtk_file = [filepath,'/Undulated_',num2str(fiber_num),'_',deformation_type];
    
    % Generate Mesh
    [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Undulated_Network_Mesh_Generator_3D( Segment,Segment_Period,Junction,Rotation,ele_size,order,period,Amp,trim_factor,'of' );

end



%% Material Properties.
step=1;
r=0.001;
Area=pi*(r)^2;E=10^6;v=0.2;G=0.5*(1+v)*E;I2=(pi/4)*r^4;I3=I2;Ip=I3+I2;% C/S area,Young's Modulus & Density
Mat(1)=Area;Mat(2)=E;Mat(3)=G;Mat(4)=I2;Mat(5)=I3;Mat(6)=Ip;
P=CPTS(1:4,:);
[XP,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,P,Patch,A0,Mat,ele_size);

% Group each increment into an individual cell for parallel loop
num_increments = size(CPTS,1)/4;
CPTS_cell = cell(num_increments);
for i=1:num_increments
    CPTS_cell{i}=CPTS(4*(i-1)+1:4*(i-1)+4,:);
end

% Find first arc - we exclude these arcs from the energy calculations
first_arc = find_first_arc(Patch);

%% Total energy

% XP gives the coordinates of control points in their reference config
% XQ gives the coordaintes of control points in their current config
% In each loop, we update the value of XQ to get the coordinates in
% the ith increment, and use those to calculate the displacements
% off the fibers in Write_VTK_Network_Line.

ME_network=zeros(size(CPTS,1)/4,1);
BE_network=ME_network;
TE_network=ME_network;
SE_network=ME_network;

% Option to use parallel loop (parfor) if desired - uncomment below lines
% parpool('local',8);
% par
for i=1:size(CPTS,1)/4%i=1:size(CPTS,1)/4
    Q=CPTS(4*(i-1)+1:4*(i-1)+4,:);
    [XQ,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,Q,Patch,A0,Mat,ele_size,first_arc);
    t1=0;t2=0;t3=0;t4=0;
    for j=1:size(ME,2)
        t1=t1+sum(ME{j});t2=t2+sum(BE{j});
        t3=t3+sum(TE{j});t4=t4+sum(SE{j});
    end
    ME_network(i,1)=t1;BE_network(i,1)=t2;TE_network(i,1)=t3;SE_network(i,1)=t4;
end
% p = gcp;
% delete(gcp('nocreate'))

%% Get strain

if type == 1
    dataname = ['forcedisp\network_',num2str(fiber_num),'_',num2str(Case),'_',deformation_type];
    savename = ['energy\network_',num2str(fiber_num),'_',num2str(Case),'_',deformation_type];
else
    dataname = ['forcedisp\network_',num2str(fiber_num),'_',num2str(Case),'_',num2str(New_Case),'_',deformation_type];
    savename = ['energy\network_',num2str(fiber_num),'_',num2str(Case),'_',num2str(New_Case),'_',deformation_type];
end

load(dataname)
strain = displacement*100;

% Plot total energies
figure; hold on;
plot(strain,ME_network)
plot(strain,BE_network)
plot(strain,TE_network)
plot(strain,SE_network)
legend('Membrane','Bending','Torsion','Total')

% Plot energy ratios
figure; hold on;
plot(strain(2:end),ME_network(2:end)./SE_network(2:end))
plot(strain(2:end),BE_network(2:end)./SE_network(2:end))
plot(strain(2:end),TE_network(2:end)./SE_network(2:end))
xlim([strain(2) strain(end)])
legend('Membrane','Bending','Torsion')
save(savename,'ME_network','BE_network','TE_network','SE_network')

end


%% Function to find when index of first arc to exclude from calculations

function arc_start = find_first_arc(Patch)
    arc_start = 0;

    for ii = 1:(size(Patch,2)-1)
        cpt = Patch(ii).cpt;
        cpt_next = Patch(ii+1).cpt;
        if ~arc_start 
            if size(cpt,2) == 4 && size(cpt_next,2) == 4
                arc_start = ii;
            end
        else
            if size(cpt,2) ~= 4
                arc_start = 0;
            end
        end            
    end
    
end