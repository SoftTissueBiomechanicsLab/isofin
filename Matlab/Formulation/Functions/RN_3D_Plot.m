function [] = RN_3D_Plot(Segment)
%% Code
figure
for i=1:size(Segment,1)
    plot3([Segment(i,1),Segment(i,4)],[Segment(i,2),Segment(i,5)],[Segment(i,3),Segment(i,6)],'g','LineWidth',2);hold on
    for j=1:2
        [value,FaceOut] = PtInBox(Segment(i,3*(j-1)+1:3*(j-1)+3));
        if value==0
% %             disp('Point on the Box!')
% %             Segment(i,3*(j-1)+1:3*(j-1)+3)
            plot3(Segment(i,3*(j-1)+1),Segment(i,3*(j-1)+2),Segment(i,3*(j-1)+3),'.b','MarkerSize',20);hold on
        else
            plot3(Segment(i,3*(j-1)+1),Segment(i,3*(j-1)+2),Segment(i,3*(j-1)+3),'.k','MarkerSize',20);hold on
        end
    end
end
hold on;
% Bottom horizontal edges
plot3([1,1],[0,1],[0,0],'k','LineWidth',1);hold on;
plot3([0,0],[0,1],[0,0],'k','LineWidth',1);hold on;
plot3([0,1],[1,1],[0,0],'k','LineWidth',1);hold on;
plot3([0,1],[0,0],[0,0],'k','LineWidth',1);hold on;
% Top horizontal edges
plot3([1,0],[1,1],[1,1],'k','LineWidth',1);hold on;
plot3([0,1],[0,0],[1,1],'k','LineWidth',1);hold on;
plot3([1,1],[0,1],[1,1],'k','LineWidth',1);hold on;
plot3([0,0],[1,0],[1,1],'k','LineWidth',1);hold on;
% Vetrical edges
plot3([0,0],[0,0],[0,1],'k','LineWidth',1);hold on;
plot3([1,1],[0,0],[0,1],'k','LineWidth',1);hold on;
plot3([1,1],[1,1],[0,1],'k','LineWidth',1);hold on;
plot3([0,0],[1,1],[0,1],'k','LineWidth',1);hold on;
xlabel('X Axis');ylabel('Y Axis');zlabel('Z Axis');
end

