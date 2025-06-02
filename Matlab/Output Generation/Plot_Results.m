%% Plot force vs displacment and save for further analysis
clearvars; close all; clc;
warning off

%% Add paths

addpath('..\Formulation\German Formulation\C_files');
addpath('..\Formulation\German Formulation\nurbs_toolbox');
addpath('..\Formulation\\German Formulation\With Torsion\Functions_RN Parallel\3D_New');
addpath('..\Input File Generation\Functions');

%% Pick fiber number, and use matching beam order and power
type = 1;% 1 for 3D straight and 2 for 3D undulated and 3 for 2D straight
fiber_num = 400;

for Case = 1
New_Case=1;

%% Mesh File parameters
order = 5;
power = 4;
ele_size=0.4/(2^power);
MS=strcat('_MS',num2str(power));
end_cpt_factor=1/5;
trim_factor=0.01;

% deformation_types are SSX for simple shear, UAX for uniaxial, 
% BIAX for biaxial, and PS for pure shear representative cube.
deformation_type = 'PS';

switch type
    case 1        
        % Load network geometry
        load(['..\Networks\Data\3D\Straight\',num2str(fiber_num),'\Case',num2str(Case)]); %#ok<*UNRCH>        
        % Input File parameters
        MeshFileName=['Straight_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS];
        % Generate Mesh
        [Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs] = Network_Mesh_Generator_3D( Segment,Junction,ele_size,order,end_cpt_factor,trim_factor,'of' );
    case 2    
        % Undulation amplitude for undulated networks
        Amp=0.025;
        % Load network geometry
        load(['..\Networks\Data\3D\Undulated\',num2str(fiber_num),'\Case',num2str(Case),'_',num2str(New_Case)]);        
        % Input File parameters
        MeshFileName=['Undulated_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_',num2str(New_Case),'_P',num2str(order),MS];
        % Generate Mesh
        [Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs] = Undulated_Network_Mesh_Generator_3D( Segment,Segment_Period,Junction,Rotation,ele_size,order,period,Amp,trim_factor,'of' );  
    case 3    
        % Mesh File parameters
        % Load network geometry
        load(['..\Networks\Data\2D\Straight\',num2str(fiber_num),'\Case',num2str(Case)]);        
        % Input File parameters
        MeshFileName= ['Straight_2D_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS];
        % Generate Mesh
        [Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs] = Network_Mesh_Generator_2D( vx,vy,ele_size,order,end_cpt_factor,trim_factor,'of' );
end

% Names of results files/paths
path_c = '..\..\Results\Outputs\';
ResultFileName=strcat(MeshFileName,'_Result.txt');
ResultFilePath_c=strcat(path_c,ResultFileName);
CPTS_c=importdata(ResultFilePath_c);

ResultFileName=strcat(MeshFileName,'_Result_Reactions.txt');
ResultFilePath_c=strcat(path_c,ResultFileName);
REACT_c=importdata(ResultFilePath_c);

P=CPTS_c(1:4,:);
for i=1:size(P,2)
    if P(1,i)==1
       CPT_num = i;
    end
end
key=1;
force=[];
displacement=[];
force_y = [];
for i=1:size(CPTS_c,1)/4
    Q=CPTS_c(1+4*(i-1):4+4*(i-1),:);
    force(key)=-sum(REACT_c(4*(i-1)+1,:));
    force_y(key)=-sum(REACT_c(4*(i-1)+2,:));
    displacement(key)=abs(P(1,CPT_num) - Q(1,CPT_num));
    key=key+1;
end
if type == 1
    savename = ['forcedisp/network_',num2str(fiber_num),'_',num2str(Case),'_',deformation_type];
else
    savename = ['forcedisp/network_',num2str(fiber_num),'_',num2str(Case),'_',num2str(New_Case),'_',deformation_type];
end
save(savename,'displacement','force','force_y')

figure
hold on
plot(displacement,force,':s','LineWidth',2,'MarkerSize',5);
xlabel('Displacement')
ylabel('Force')

plot(displacement,force_y,':s','LineWidth',2,'MarkerSize',5);
legend('X','Y')
% ylim([0 ceil(max([force(end) force_y(end)]))])
end
legend




