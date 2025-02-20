function [ F_DOF, D_DOF ] = Apply_BC_Network( Cpt,Edge_D,Edge_F )
%% Code
Right=[];Left=[];Top=[];Bottom=[];
keyR=1;keyL=1;keyT=1;keyB=1;
%% Identify cpts on the edges
for i=1:size(Cpt,2)
    % bottom edge
    if Cpt(2,i)==0
       Bottom(1,keyB)=i;keyB=keyB+1;
    end
    % right edge
    if Cpt(1,i)==1
       Right(1,keyR)=i;keyR=keyR+1;
    end
    % top edge
    if Cpt(2,i)==1
        Top(1,keyT)=i;keyT=keyT+1;
    end
    % left edge
    if Cpt(1,i)==0
        Left(1,keyL)=i;keyL=keyL+1;
    end
end
%% Form D_DOF (for applying displacement boundary conditions on the edges)
D_DOF=Inf*ones(size(Cpt,2),size(Cpt,1));
% Bottom Edge
for i=1:size(Bottom,2)
    for j=1:size(Edge_D,2)
        if Edge_D(1,j)~=Inf
           D_DOF(Bottom(i),j)=Edge_D(1,j);
        end
    end
end
% Right Edge
for i=1:size(Right,2)
    for j=1:size(Edge_D,2)
        if Edge_D(2,j)~=Inf
           D_DOF(Right(i),j)=Edge_D(2,j);
        end
    end
end
% Top Edge
for i=1:size(Top,2)
    for j=1:size(Edge_D,2)
        if Edge_D(3,j)~=Inf
           D_DOF(Top(i),j)=Edge_D(3,j);
        end
    end
end
% Left Edge 
for i=1:size(Left,2)
    for j=1:size(Edge_D,2)
        if Edge_D(4,j)~=Inf
           D_DOF(Left(i),j)=Edge_D(4,j);
        end
    end
end       
%% Form F_DOF (for applying force boundary conditions on the edges)
F_DOF=zeros(size(Cpt,2),size(Cpt,1));
% Bottom Edge
for i=1:size(Bottom,2)
    for j=1:size(Edge_F,2)
           F_DOF(Bottom(i),j)=Edge_F(1,j)/size(Bottom,2);
    end
end
% Right Edge
for i=1:size(Right,2)
    for j=1:size(Edge_F,2)
           F_DOF(Right(i),j)=Edge_F(2,j)/size(Right,2);
    end
end
% Top Edge
for i=1:size(Top,2)
    for j=1:size(Edge_F,2)
           F_DOF(Top(i),j)=Edge_F(3,j)/size(Top,2);
    end
end
% Left Edge 
for i=1:size(Left,2)
    for j=1:size(Edge_F,2)
           F_DOF(Left(i),j)=Edge_F(4,j)/size(Left,2);
    end
end
end

