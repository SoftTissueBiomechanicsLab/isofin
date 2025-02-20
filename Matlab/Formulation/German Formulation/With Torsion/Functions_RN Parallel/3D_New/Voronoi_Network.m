function [ vx,vy ] = Voronoi_Network( npt,fig )
%Voronoi_Network : Creates Voronoi network in a unit square.
%% Code
X=rand(npt,2);
Bound=[0,0;1,0;1,1;0,1];
[V,C,XY]=VoronoiLimit(X(:,1),X(:,2),'bs_ext',Bound,'figure','off');
V=round(V,6);
vx=[];vy=[];key=1;
%% Find interior vetrices
for i=1:size(C,1)
    c=C{i};
    for j=1:size(c,1)-1
       x1= V(c(j),1);x2=V(c(j+1),1);y1 = V(c(j),2);y2 = V(c(j+1),2);
       if (y1==0&&y2==0)||(y1==1&&y2==1)||(x1==0&&x2==0)||(x1==1&&x2==1) % eliminate vertices on the boundary
       else
           vx(1,key) = x1;vy(1,key) = y1;
           vx(2,key) = x2;vy(2,key) = y2;
           key=key+1;
       end
    end
   x1= V(c(j+1),1);x2=V(c(1),1);y1 = V(c(j+1),2);y2 = V(c(1),2);
   if (y1==0&&y2==0)||(y1==1&&y2==1)||(x1==0&&x2==0)||(x1==1&&x2==1) % eliminate vertices on the boundary
   else
       vx(1,key) = x1;vy(1,key) = y1;
       vx(2,key) = x2;vy(2,key) = y2;
       key=key+1;
   end   
end
% % figure
% % plot(vx,vy,'b')
remove=[];key=1;
for i=1:size(vx,2)
    test_x=[vx(2,i);vx(1,i)];
    for j=1:size(vx,2)
        if (test_x(1)==vx(1,j) && test_x(2)==vx(2,j)) && j>i
            remove(key)=j;key=key+1;
        end
    end
end
vx(:,remove)=[];vy(:,remove)=[];
if fig=='on'
    figure
    plot(vx,vy,'r')
end
end

