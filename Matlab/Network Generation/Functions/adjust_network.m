function [Network,junction_number] = adjust_network(Network)

% Find how many fibers are coming out each Junction
junction_number = how_many_fibers(Network);    
% Remove junctions that have 1 or 0 fibers
zero_ind = junction_number == 0;
one_ind = junction_number == 1;
Network.Junction(zero_ind|one_ind,:) = [];
% Get new junction_number for debugging purposes if desired
junction_number = how_many_fibers(Network);