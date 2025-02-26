 %% TXT file writing and reading.
clearvars;close all;clc
tic
%% Add paths
addpath('../Formulation/German Formulation/C_files')
addpath('../Formulation/German Formulation/nurbs_toolbox')
addpath('../Formulation/German Formulation/With Torsion/Functions_RN Parallel/3D_New');
addpath('Functions')
%% 3D Network :
type=2; % 1 for 3D straight and 2 for 3D undulated and 3 for 2D straight
fiber_num=350;
Case=1;
New_Case=1;
% deformation_types are SSX for simple shear, UAX for uniaxial, and BIAX for biaxial.
% You can adjust these as desired.
deformation_type = 'BIAX';

for power=4
    if type==1
        %% Mesh File parameters
        ele_size=0.4/(2^power);
        MS=strcat('_MS',num2str(power));
        order=5;
        end_cpt_factor=1/5;
        trim_factor=0.01;
        % Load network geometry        
        filename = ['../Networks/Data/3D/Straight/',num2str(fiber_num),'/Case',num2str(Case)];
        load(filename)

        %% Input File parameters
        MeshFileName = ['Straight_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS];
        InputFileName = ['input_S_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS];
    elseif type==2
        %% Mesh File parameters
        ele_size=0.4/(2^power);
        MS=strcat('_MS',num2str(power));
        order=5;
        trim_factor=0.01;
        Amp=0.05;
        % Load network geometry
        filename = ['../Networks/Data/3D/Undulated/',num2str(fiber_num),'/Case',num2str(Case),'_',num2str(New_Case)];
        load(filename)

        %% Input File parameters
        MeshFileName=strcat('Undulated_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_',num2str(New_Case),'_P',num2str(order),MS);
        InputFileName=strcat('input_U_',deformation_type,'_',num2str(fiber_num),'_Case',num2str(Case),'_',num2str(New_Case),'_P',num2str(order),MS);
    elseif type==3
        %% Mesh File parameters
        ele_size=0.4/(2^power);MS=strcat('_MS',num2str(power));
        order=5;end_cpt_factor=0.2;trim_factor=0.02;
        % load network geometry
        filename = ['../Networks/Data/2D/Straight/',num2str(fiber_num),'/Case',num2str(Case)];
        load(filename)        

        %% Input File parameters
        MeshFileName=strcat('Straight_2D_BAT_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS);
        % InputFileName=strcat(MeshFileName,'_Input');
        InputFileName=strcat('input_S_2D_',num2str(fiber_num),'_Case',num2str(Case),'_P',num2str(order),MS);
    end
    % Material Properties
    r=0.001;E=10^6;v=0.2;
    % Boundary conditions
    Max_Disp=0.1;
    % Big number to identify unassigned displacement dofs.
    u0=10^10;
    % Displacement BC - can switch what faces BCs are applied to as needed
    switch deformation_type
        case 'SSX'
        Face_D(1,:)=[u0,u0,u0,u0]; % Face 1 (x=1)
        Face_D(2,:)=[u0,u0,u0,u0]; % Face 2 (x=0)
        Face_D(3,:)=[u0,u0,u0,u0]; % Face 3 (y=1)
        Face_D(4,:)=[u0,u0,u0,u0]; % Face 4 (y=0)
        Face_D(5,:)=[Max_Disp,0,0,0]; % Face 5 (z=1)
        Face_D(6,:)=[0,0,0,0]; % Face 6 (z=0)
        case 'UAX'
        Face_D(1,:)=[u0,u0,u0,u0]; % Face 1 (x=1)
        Face_D(2,:)=[u0,u0,u0,u0]; % Face 2 (x=0)
        Face_D(3,:)=[u0,u0,u0,u0]; % Face 3 (y=1)
        Face_D(4,:)=[u0,u0,u0,u0]; % Face 4 (y=0)
        Face_D(5,:)=[0,0,Max_Disp,0]; % Face 5 (z=1)
        Face_D(6,:)=[0,0,0,0]; % Face 6 (z=0)
        case 'BIAX'
        Face_D(1,:)=[Max_Disp,0,0,0]; % Face 1 (x=1)
        Face_D(2,:)=[0,0,0,0]; % Face 2 (x=0)
        Face_D(3,:)=[u0,u0,u0,u0]; % Face 3 (y=1)
        Face_D(4,:)=[u0,u0,u0,u0]; % Face 4 (y=0)
        Face_D(5,:)=[0,0,Max_Disp,0]; % Face 5 (z=1)
        Face_D(6,:)=[0,0,0,0]; % Face 6 (z=0)

    end
    % Force BC
    Face_F(1,:)=[0,0,0,0];% Face 1 (x=1)
    Face_F(2,:)=[0,0,0,0]; % Face 2 (x=0)
    Face_F(3,:)=[0,0,0,0]; % Face 3 (y=1)
    Face_F(4,:)=[0,0,0,0]; % Face 4 (y=0)
    Face_F(5,:)=[0,0,0,0]; % Face 5 (z=1)
    Face_F(6,:)=[0,0,0,0]; % Face 6 (z=0)
    % Analysis parameters
    num_threads=36; % No of cores requested for analysis
    num_inc=200; % No of increments.
    max_numiter=50; % Max no of iterations per increment
    max_attempts=10; % Max attempts per increment. (applied displacement halved per attempt)
    tol_D=10^-2;tol_R=5*10^-3; % tolerence for displacement and residue criterion.
    
    %% plot network
    % RN_3D_Plot(Segment)
    
    %% Generate Mesh
    if type==1
        [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator_3D( Segment,Junction,ele_size,order,end_cpt_factor,trim_factor,'of' );
    end
    if type==2
        [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Undulated_Network_Mesh_Generator_3D( Segment,Segment_Period,Junction,Rotation,ele_size,order,period,Amp,trim_factor,'of' );
    end
    if type==3
        [ Cpt,Wts,Patch,Patch_Bs,A0,A0_Bs ] = Network_Mesh_Generator_2D( vx,vy,ele_size,order,end_cpt_factor,trim_factor,'of' );
    end


%% Write input file data to Text file
path = '../../Input_files/';
% filename=strcat(path,'Trial_3D','.txt');
filename=strcat(path,InputFileName,'.txt');
fileID = fopen(filename,'w');
% Mesh file name
fprintf(fileID,strcat(MeshFileName,'\r\n'));
% Material Properties
fprintf(fileID,'%f %f %f \r\n',r,E,v);
% Boundary Conditions
fprintf(fileID,'%f \r\n',Max_Disp);
fprintf(fileID,'%f \r\n',u0);
for i=1:size(Face_D,1)
fprintf(fileID,'%f %f %f %f \r\n',Face_D(i,:));
end
for i=1:size(Face_F,1)
fprintf(fileID,'%f %f %f %f \r\n',Face_F(i,:));
end
% Analysis Parameters
fprintf(fileID,'%d %d %d %d \r\n',num_threads,num_inc,max_numiter,max_attempts);
fprintf(fileID,'%f %f \r\n',tol_D,tol_R);
fclose(fileID);


%% Write mesh file data to Text file

path = '../../Mesh_files/';
% filename=strcat(path,'Trial_3D','.txt');
filename=strcat(path,MeshFileName,'.txt');
fileID = fopen(filename,'w');
num_cpt=size(Cpt,2);
% Control Points
fprintf(fileID,'%d \r\n',num_cpt);
for i=1:size(Cpt,2)
fprintf(fileID,'%f %f %f %f \r\n',Cpt(:,i));
end
% Weights
fprintf(fileID,'%f ',Wts);fprintf(fileID,'\r\n');
% Patch
fprintf(fileID,'%d ',size(Patch,2));fprintf(fileID,'\r\n');
for i = 1:size(Patch,2)
    fprintf(fileID,'%d %d %d ',size(Patch(i).kv,2),size(Patch(i).wts,2),Patch(i).order);
    fprintf(fileID,'\r\n');
    fprintf(fileID,'%f ',Patch(i).kv);fprintf(fileID,'\r\n');
    fprintf(fileID,'%f ',Patch(i).wts);fprintf(fileID,'\r\n');
    fprintf(fileID,'%d ',Patch(i).num_cpt);fprintf(fileID,'\r\n');
    for j=1:size(Patch(i).cpt,2)
        fprintf(fileID,'%f %f %f %f \r\n',Patch(i).cpt(:,j));
    end
end
% Patch_Bs
fprintf(fileID,'%d ',size(Patch_Bs,2));fprintf(fileID,'\r\n');
for i = 1:size(Patch_Bs,2)
    fprintf(fileID,'%d %d %d ',size(Patch_Bs(i).kv,2),size(Patch_Bs(i).wts,2),Patch_Bs(i).order);
    fprintf(fileID,'\r\n');
    fprintf(fileID,'%f ',Patch_Bs(i).kv);fprintf(fileID,'\r\n');
    fprintf(fileID,'%f ',Patch_Bs(i).wts);fprintf(fileID,'\r\n');
    fprintf(fileID,'%d ',Patch_Bs(i).num_cpt);fprintf(fileID,'\r\n');
    for j=1:size(Patch_Bs(i).cpt,2)
        fprintf(fileID,'%f %f %f %f \r\n',Patch_Bs(i).cpt(:,j));
    end
end
% A0
for i=1:size(Patch,2)
fprintf(fileID,'%f %f %f \r\n',A0(1,3*(i-1)+1:3*(i-1)+3));
fprintf(fileID,'%f %f %f \r\n',A0(2,3*(i-1)+1:3*(i-1)+3));
fprintf(fileID,'%f %f %f \r\n',A0(3,3*(i-1)+1:3*(i-1)+3));
end
% A0_Bs
for i=1:size(Patch_Bs,2)
fprintf(fileID,'%f %f %f \r\n',A0_Bs(1,3*(i-1)+1:3*(i-1)+3));
fprintf(fileID,'%f %f %f \r\n',A0_Bs(2,3*(i-1)+1:3*(i-1)+3));
fprintf(fileID,'%f %f %f \r\n',A0_Bs(3,3*(i-1)+1:3*(i-1)+3));
end
fclose(fileID);
end







