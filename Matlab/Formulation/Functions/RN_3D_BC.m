function [ F_DOF, D_DOF ] = RN_3D_BC( Cpt,Face_D,Face_F )
%% Apply boundary conditions on a unit cube.
%% Code
Face_x1=[];Face_x0=[];Face_y1=[];Face_y0=[];Face_z1=[];Face_z0=[]; % Stores ID of cpts on these faces. (Face_x1  = Face perpendicular to x at x=1)
key=ones(1,6);

%% Identify cpts on the faces
for i=1:size(Cpt,2)
    % face 1 (x=1)
    if Cpt(1,i)==1 && sum(ismember(Cpt(2,i),[0,1]))==0 && sum(ismember(Cpt(3,i),[0,1]))==0
        Face_x1(1,key(1))=i;key(1)=key(1)+1;
    end
    % face 2 (x=0)
    if Cpt(1,i)==0 && sum(ismember(Cpt(2,i),[0,1]))==0 && sum(ismember(Cpt(3,i),[0,1]))==0
        Face_x0(1,key(2))=i;key(2)=key(2)+1;
    end
    % face 3 (y=1)
    if Cpt(2,i)==1 && sum(ismember(Cpt(1,i),[0,1]))==0 && sum(ismember(Cpt(3,i),[0,1]))==0
        Face_y1(1,key(3))=i;key(3)=key(3)+1;
    end
    % face 4 (y=0)
    if Cpt(2,i)==0 && sum(ismember(Cpt(1,i),[0,1]))==0 && sum(ismember(Cpt(3,i),[0,1]))==0
        Face_y0(1,key(4))=i;key(4)=key(4)+1;
    end
    % face 5 (z=1)
    if Cpt(3,i)==1 && sum(ismember(Cpt(1,i),[0,1]))==0 && sum(ismember(Cpt(2,i),[0,1]))==0
        Face_z1(1,key(5))=i;key(5)=key(5)+1;
    end
    % face 6 (z=0)
    if Cpt(3,i)==0 && sum(ismember(Cpt(1,i),[0,1]))==0 && sum(ismember(Cpt(2,i),[0,1]))==0
        Face_z0(1,key(6))=i;key(6)=key(6)+1;
    end
    
end

%% Form D_DOF (for applying displacement boundary conditions on the edges)
D_DOF=Inf*ones(size(Cpt,2),size(Cpt,1));
% face 1 (x=1)
for i=1:size(Face_x1,2)
    for j=1:size(Face_D,2)
        if Face_D(1,j)~=Inf
           D_DOF(Face_x1(i),j)=Face_D(1,j);
        end
    end
end
% face 2 (x=0)
for i=1:size(Face_x0,2)
    for j=1:size(Face_D,2)
        if Face_D(2,j)~=Inf
           D_DOF(Face_x0(i),j)=Face_D(2,j);
        end
    end
end  
% face 3 (y=1)
for i=1:size(Face_y1,2)
    for j=1:size(Face_D,2)
        if Face_D(3,j)~=Inf
           D_DOF(Face_y1(i),j)=Face_D(3,j);
        end
    end
end
% face 4 (y=0)
for i=1:size(Face_y0,2)
    for j=1:size(Face_D,2)
        if Face_D(4,j)~=Inf
           D_DOF(Face_y0(i),j)=Face_D(4,j);
        end
    end
end
% face 5 (z=1)
for i=1:size(Face_z1,2)
    for j=1:size(Face_D,2)
        if Face_D(5,j)~=Inf
           D_DOF(Face_z1(i),j)=Face_D(5,j);
        end
    end
end
% face 6 (z=0)
for i=1:size(Face_z0,2)
    for j=1:size(Face_D,2)
        if Face_D(6,j)~=Inf
           D_DOF(Face_z0(i),j)=Face_D(6,j);
        end
    end
end

%% Form F_DOF (for applying force boundary conditions on the edges)
F_DOF=zeros(size(Cpt,2),size(Cpt,1));
% face 1 (x=1)
for i=1:size(Face_x1,2)
    for j=1:size(Face_F,2)
%         if Face_F(1,j)~=Inf
           F_DOF(Face_x1(i),j)=Face_F(1,j)/size(Face_x1,2);
%         end
    end
end
% face 2 (x=0)
for i=1:size(Face_x0,2)
    for j=1:size(Face_F,2)
%         if Face_F(2,j)~=Inf
           F_DOF(Face_x0(i),j)=Face_F(2,j)/size(Face_x0,2);
%         end
    end
end  
% face 3 (y=1)
for i=1:size(Face_y1,2)
    for j=1:size(Face_F,2)
%         if Face_F(3,j)~=Inf
           F_DOF(Face_y1(i),j)=Face_F(3,j)/size(Face_y1,2);
%         end
    end
end
% face 4 (y=0)
for i=1:size(Face_y0,2)
    for j=1:size(Face_F,2)
%         if Face_F(4,j)~=Inf
           F_DOF(Face_y0(i),j)=Face_F(4,j)/size(Face_y0,2);
%         end
    end
end
% face 5 (z=1)
for i=1:size(Face_z1,2)
    for j=1:size(Face_F,2)
%         if Face_F(5,j)~=Inf
           F_DOF(Face_z1(i),j)=Face_F(5,j)/size(Face_z1,2);
%         end
    end
end
% face 6 (z=0)
for i=1:size(Face_z0,2)
    for j=1:size(Face_F,2)
%         if Face_F(6,j)~=Inf
           F_DOF(Face_z0(i),j)=Face_F(6,j)/size(Face_z0,2);
%         end
    end
end
end

