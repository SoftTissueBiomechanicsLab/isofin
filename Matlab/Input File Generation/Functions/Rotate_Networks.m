function [Segment, Junction] = Rotate_Networks(Segment,Junction,rotation,side_length)
% Use this function to rotate your

Segment_1 = zeros(size(Segment,1),3);
Segment_2 = zeros(size(Segment,1),3);
% Only first three columns of junction are used
New_Junction = zeros(size(Junction,1),3);

switch rotation
    case 1
        % Do nothing
    case 2
        % Rotate 90 degrees around z-axis (Right face)
        R = [0 -1 0
             1 0 0
             0 0 1];

            for i = 1:size(Segment,1)
                Segment_1(i,:) = Segment(i,1:3)*R;
                Segment_2(i,:) = Segment(i,4:6)*R;   
            end
            for i = 1:size(Junction,1)
                New_Junction(i,:) = Junction(i,1:3)*R;
            end
        Segment = [Segment_1 Segment_2];
        Junction = New_Junction;

    case 3
        % Rotate 180 degrees around z-axis (Back face)
         R = [0 -1 0
             1 0 0
             0 0 1];

            for i = 1:size(Segment,1)
                Segment_1(i,:) = Segment(i,1:3)*R*R;
                Segment_2(i,:) = Segment(i,4:6)*R*R;   
            end
            for i = 1:size(Junction,1)
                New_Junction(i,:) = Junction(i,1:3)*R;
            end
        Segment = [Segment_1 Segment_2];
        Junction = New_Junction;


    case 4
        % Rotate 270 degrees around z-axis (Left face)
         R = [0 -1 0
             1 0 0
             0 0 1];

            for i = 1:size(Segment,1)
                Segment_1(i,:) = Segment(i,1:3)*R*R*R;
                Segment_2(i,:) = Segment(i,4:6)*R*R*R;   
            end
            for i = 1:size(Junction,1)
                New_Junction(i,:) = Junction(i,1:3)*R;
            end
        Segment = [Segment_1 Segment_2];
        Junction = New_Junction;

    case 5
        % Rotate 90 degress around x-axis (Top face)
        R = [1 0 0
             0 0 -1
             0 1 0];
        for i = 1:size(Segment,1)
            Segment_1(i,:) = Segment(i,1:3)*R;
            Segment_2(i,:) = Segment(i,4:6)*R;   
        end
        for i = 1:size(Junction,1)
                New_Junction(i,:) = Junction(i,1:3)*R;
        end
        Segment = [Segment_1 Segment_2];
        Junction = New_Junction;

    case 6
        % Rotate 270 degress around x-axis (Bottom face)
        R = [1 0 0
             0 0 1
             0 -1 0];
        for i = 1:size(Segment,1)
            Segment_1(i,:) = Segment(i,1:3)*R;
            Segment_2(i,:) = Segment(i,4:6)*R;   
        end
        for i = 1:size(Junction,1)
                New_Junction(i,:) = Junction(i,1:3)*R;
        end
        Segment = [Segment_1 Segment_2];
        Junction = New_Junction;

end

% Moves network back to original 
min_x = min(Segment(:,[1 4]),[],'all');
min_y = min(Segment(:,[2 5]),[],'all');
min_z = min(Segment(:,[3 6]),[],'all');

if min_x < 0
    Segment(:,[1 4]) = Segment(:,[1 4]) + side_length;
    Junction(:,1) = Junction(:,1) + side_length;
end
if min_y < 0
    Segment(:,[2 5]) = Segment(:,[2 5]) + side_length;
    Junction(:,2) = Junction(:,2) + side_length;

end
if min_z < 0
    Segment(:,[3 6]) = Segment(:,[3 6]) + side_length;
    Junction(:,3) = Junction(:,3) + side_length;
end
