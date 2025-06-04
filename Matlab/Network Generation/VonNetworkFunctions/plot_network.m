function []=plot_network(Segment,color_l,fig_num)
lw = 0.2;
n_fibers = size(Segment,1);

figure(fig_num)
hold on
for i = 1:n_fibers

    start_point = Segment(i,1:3);
    end_point = Segment(i,4:6);

       XYZ = [start_point;end_point];

    plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'-','color',color_l,'linewidth',lw)

end
axis equal
view([-36 27])
grid on
hold off