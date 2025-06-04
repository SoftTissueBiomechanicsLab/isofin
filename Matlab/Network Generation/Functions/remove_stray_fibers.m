function Segment = remove_stray_fibers(Segment)

%% Stray fibers are fibers that are not connected to another fiber

% Loop through Segment. Create row ids to loop through. If a fiber node
% appears in multiple rows, remove the row_id. If not, keep it in and
% remove the row after.

row_id = 1:1:size(Segment,1);
keep_id = [];

for ii = row_id
    node_one = Segment(ii,1:3);
    node_two = Segment(ii,4:6);

    % group_one is all nodes in columns 1:3 excluding the current fiber
    % being compared to these rows
    group_one = Segment;
    group_one(:,4:6) = []; 
    group_one(row_id(ii),:) = [];

    % group_two is all nodes in columns 4:6 excluding the current fiber
    % being compared to these rows
    group_two = Segment;
    group_two(:,1:3) = []; 
    group_two(row_id(ii),:) = [];

    % col_ind(row#,col#) is an index to see if any other nodes in group_col# are
    % the same as fiber_row#
    col_ind = cell(2,2);
    col_ind{1,1} = node_one == group_one;
    col_ind{1,2} = node_one == group_two;
    col_ind{2,1} = node_two == group_one;
    col_ind{2,2} = node_two == group_two;

    for j = 1:numel(col_ind)
        temp = col_ind{j};
        node_check = sum(temp,2);
        if sum(node_check == 3)
            keep_id = [keep_id ii];
            break
        end
    end
end

Segment = Segment(keep_id,:);