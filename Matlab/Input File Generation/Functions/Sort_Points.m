function [sorted_points,sort_id] = Sort_Points(points)
% points is an array holding x,y,z coordinates of a node in each column

% center = mean(points);
% theta = cart2pol(points(:,1)-center(1),points(:,2)-center(2));
% theta_wrapped = mod(theta,2*pi);
% theta_wrapped(theta_wrapped==0 & theta>0) = 2*pi;
% [~, sort_id] = sort(theta_wrapped, 'descend');
% % Swap values in sort_id - works for four and five fiber junctions 
% id_one = sort_id==1;
% id_two = sort_id==2;
% id_three = sort_id==3;
% id_four = sort_id==4;
% sort_id(id_two) = 1;
% sort_id(id_one) = 2;
% sort_id(id_four) = 3;
% sort_id(id_three) = 4;
% sorted_points = points(sort_id,:);

% Find the bottom-left point (min y, then min x if tie)
[~, ref_id] = min(points(:,2) + 1e-6 * points(:,1));  % break ties using x
refPoint = points(ref_id,:);

% Compute angles from the reference point to all others
vectors = points - refPoint;  % Shift to origin at refPoint
angles = atan2(vectors(:,2), vectors(:,1));  % Compute angle in radians

% Normalize angles to [0, 2*pi)
angles = mod(angles, 2*pi);

% Sort by angle
[~, sort_id] = sort(angles);

% Ensure the reference point is first in the order
% (optional, depending on use-case; remove this if unnecessary)
refSortPos = find(sort_id == ref_id, 1);
sort_id = circshift(sort_id, -refSortPos + 1);

% Apply sorting
sorted_points = points(sort_id, :);

end
