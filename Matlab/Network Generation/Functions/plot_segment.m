function [] = plot_segment(Segment,varargin)

if nargin > 1
    angles = varargin{1};
else
    angles = [-37.5,30];
end

figure
hold on
for i = 1:size(Segment,1)
    plot3(Segment(i,[1 4]),Segment(i,[2 5]),Segment(i,[3 6]))
end
axis equal
xlabel('X');ylabel('Y');zlabel('Z');

view(angles(1),angles(2))
