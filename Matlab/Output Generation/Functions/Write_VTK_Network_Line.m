function [] = Write_VTK_Network_Line(XP,XQ,EPS,KAPPA,ME,BE,TE,SE,filepath,step)
 

if step < 10
    extension = ['000',num2str(step)];
elseif step < 100 && step > 9
    extension = ['00',num2str(step)];
else
    extension = ['0',num2str(step)];
end
 
max=size(XP,2);
npoints=0;
for j=1:max 
    npoints=npoints+size(XP{j},1);
end

fid1 = fopen(strcat(filepath,extension,'.vtk'),'w');
 

fprintf(fid1,'%s\n','# vtk DataFile Version 3.0');
fprintf(fid1,'%s\n','test_file');
fprintf(fid1,'%s\n','ASCII');
fprintf(fid1,'%s\n','DATASET POLYDATA');
%% Points
fprintf(fid1,'POINTS %i float\n',npoints);
num_lines=0;Point_ID={};key=0;
for j=1:max 
    LQ=XQ{j};pt_ID=[];
    for i = 1 : size(LQ,1)
%         fprintf(fid1,'%10.5f %10.5f %10.5f \n',[LQ(i,1),LQ(i,2),LQ(i,3)]);
        fprintf(fid1,'%f %f %f \n',[LQ(i,1),LQ(i,2),LQ(i,3)]);
        pt_ID(1,i)=key;key=key+1;
    end
    Point_ID{j}=pt_ID;
    num_lines=size(LQ,1)-1+num_lines;
end
%% Lines
fprintf(fid1,'LINES %i %i\n',num_lines, (3*num_lines));
for j=1:max
   LQ=XQ{j};pt_ID=Point_ID{j};
   for i = 1 : size(pt_ID,2)-1
        fprintf(fid1,'%i %i %i\n',2,pt_ID(i),pt_ID(i+1));
   end
end
% Color Coding
fprintf(fid1,'POINT_DATA %i\n', npoints);
% X displacement
fprintf(fid1,'SCALARS U1 float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
num_u1=0;
for j=1:size(XP,2)
    LP=XP{j};
    for i = 1 : size(LP,1)
        num_u1=num_u1+size(LP,1);
    end
end
key=1;U1=zeros(1,num_u1);
U2=U1;U3=U1;U=U1;R1=U1;E=U1;R2=U1;R3=U1;me=U1;be=U1;te=U1;se=U1;
for j=1:size(XP,2)
    LP=XP{j};LQ=XQ{j};
    for i = 1 : size(LP,1)
        U1(1,key)=LQ(i,1)-LP(i,1);key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',U1(1,i));
end
% Y displacement
fprintf(fid1,'SCALARS U2 float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
key=1;
for j=1:size(XP,2)
    LP=XP{j};LQ=XQ{j};
    for i = 1 : size(LP,1)
        U2(1,key)=LQ(i,2)-LP(i,2);key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',U2(1,i));
end

% Z displacement
fprintf(fid1,'SCALARS U3 float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
key=1;
for j=1:size(XP,2)
    LP=XP{j};LQ=XQ{j};
    for i = 1 : size(LP,1)
        U3(1,key)=LQ(i,3)-LP(i,3);key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',U3(1,i));
end

% Torsion
fprintf(fid1,'SCALARS R1 float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
key=1;
for j=1:size(XP,2)
    LP=XP{j};LQ=XQ{j};
    for i = 1 : size(LP,1)
        R1(1,key)=LQ(i,4)-LP(i,4);key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',R1(1,i));
end

% Displacement Magnitude
fprintf(fid1,'SCALARS U float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
key=1;
for j=1:size(XP,2)
    LP=XP{j};LQ=XQ{j};
    for i = 1 : size(LP,1)
        U(1,key)=norm(LQ(i,:)-LP(i,:));key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',U(1,i));
end

% Membrane Strain
fprintf(fid1,'SCALARS E float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
key=1;
for j=1:size(XP,2)
    LP=EPS{j};
    for i = 1 : size(LP,1)
        E(1,key)=LP(i);key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',E(1,i));
end

% Bending Curvatures
fprintf(fid1,'SCALARS R2 float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
key=1;
for j=1:size(XP,2)
    LP=KAPPA{j};
    for i = 1 : size(LP,1)
        R2(1,key)=LP(i,1);R3(1,key)=LP(i,2);
        key=key+1;
    end
end
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',R2(1,i));
end

fprintf(fid1,'SCALARS R3 float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',R3(1,i));
end

% Energies (ME,BE,TE)
key=1;
for j=1:size(XP,2)
    LP1=ME{j};LP2=BE{j};LP3=TE{j};LP4=SE{j};
    for i = 1 : size(LP1,1)
        me(1,key)=LP1(i,1);be(1,key)=LP2(i,1);
        te(1,key)=LP3(i,1);se(1,key)=LP4(i,1);
        key=key+1;
    end
end
me_f=me./se;be_f=be./se;te_f=te./se;

% fprintf(fid1,'SCALARS me_f float\n');
% fprintf(fid1,'LOOKUP_TABLE default\n');
% for i = 1 : npoints
%     fprintf(fid1,'%14.7f\n',me_f(1,i));
% end
% 
% fprintf(fid1,'SCALARS be_f float\n');
% fprintf(fid1,'LOOKUP_TABLE default\n');
% for i = 1 : npoints
%     fprintf(fid1,'%14.7f\n',be_f(1,i));
% end
% 
% fprintf(fid1,'SCALARS te_f float\n');
% fprintf(fid1,'LOOKUP_TABLE default\n');
% for i = 1 : npoints
%     fprintf(fid1,'%14.7f\n',te_f(1,i));
% end

fprintf(fid1,'SCALARS me float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',me(1,i));
end

fprintf(fid1,'SCALARS be float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',be(1,i));
end

fprintf(fid1,'SCALARS te float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',te(1,i));
end

fprintf(fid1,'SCALARS se float\n');
fprintf(fid1,'LOOKUP_TABLE default\n');
for i = 1 : npoints
    fprintf(fid1,'%14.7f\n',se(1,i));
end

fclose(fid1);
end