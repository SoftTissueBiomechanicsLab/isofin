function junction_number = how_many_fibers(Network)
junction_number = zeros(size(Network.Junction,1),1);
for ii = 1:length(junction_number)
    temp = Network.Junction(ii,1:3);
    check_one = sum(temp == Network.Segment(:,1:3),2);
    check_two = sum(temp == Network.Segment(:,4:6),2);
    % Gets how many fibers come out of said junction
    junction_number(ii) = sum(floor(check_one/3) + floor(check_two/3));
end