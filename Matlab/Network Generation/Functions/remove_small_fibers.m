function Segment = remove_small_fibers(Segment,varargin)
if nargin > 1
    tol = varargin{1};
else
    tol = 3e-2;
end

remove_ind = false(size(Segment,1),1);
for ii = 1:length(remove_ind)
    fiber = [Segment(ii,1)-Segment(ii,4) Segment(ii,2)-Segment(ii,5) Segment(ii,3)-Segment(ii,6)];
    fiber_length = norm(fiber);
    if fiber_length <= tol
        remove_ind(ii) = true;
    end
end
Segment(remove_ind,:) = [];