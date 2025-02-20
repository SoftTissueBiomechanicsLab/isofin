function [] = PostProcess_3D_Network(filepath)
%% Code
%% Extend dangling fibers till closest face. 
load(filepath,'Junction','Segment');
Seg_temp=Segment; Jun_temp=Junction;
for i=1:size(Seg_temp,1)
    for j=1:size(Jun_temp)
        if Seg_temp(i,1:3)==Jun_temp(j,1:3)
            key1=1;break;
        else
            key1=0;
        end
    end
    for j=1:size(Jun_temp)
        if Seg_temp(i,4:6)==Jun_temp(j,1:3)
            key2=1;break;
        else
            key2=0;
        end
    end
    if key2*key1 == 0 % Fiber is dangling!
       % edit the dangling end such that it lies on the face in the line o the fiber direction
       if key1==0 % fiber dangling from end 1
          n1=(Seg_temp(i,1:3)-Seg_temp(i,4:6))/norm(Seg_temp(i,1:3)-Seg_temp(i,4:6));
          X(:,1)=Seg_temp(i,4:6)+10*n1;
          [~,FaceOut] = PtInBox(X(:,1));
          X(:,1)=terminate(Seg_temp(i,4:6)',X(:,1),FaceOut);
          Seg_temp(i,1:3)=X(:,1);
       end
       if key2==0 % fiber dangling from end 1
           n2=(Seg_temp(i,4:6)-Seg_temp(i,1:3))/norm(Seg_temp(i,1:3)-Seg_temp(i,4:6));
           X(:,1)=Seg_temp(i,1:3)+10*n2;
           [~,FaceOut] = PtInBox(X(:,1));
           X(:,1)=terminate(Seg_temp(i,1:3)',X(:,1),FaceOut);
           Seg_temp(i,4:6)=X(:,1);               
       end
    end
end
Segment=Seg_temp;
% Adjust fiber ends close to the faces.
for i=1:size(Segment,1)
    for j=1:6
        if abs(Segment(i,j)-0) < 10^-6
           Segment(i,j)=0;
        end
        if abs(Segment(i,j)-1) < 10^-6
            Segment(i,j)=1;
        end
        
    end
end
save(filepath,'Junction','Segment')
end











