function [NSegment,FJunctions,FSegmentsCon,FJunctiosCon,FLocConNum, InnerConNum,FFibLength,theta,FOnEdge] = arrange_network(Segment,L,ref_vec)




FJunctions = unique([Segment(:,1:3);Segment(:,4:6)],'rows');
FSegmentsCon = [];

% Define Unique Segments Connectivity
FFibLength = zeros(size(Segment,1),1);
rm_IDs = [];
for i = 1:size(Segment,1)
    r0 = Segment(i,1:3);
    rf = Segment(i,4:6);
    fib = rf-r0;
    theta(i,1) = acos(dot(ref_vec,fib)/norm(fib));
    FFibLength(i,1) = norm(rf-r0);
    index=ismember(FJunctions,[r0;rf],'rows');
    IDs = find(index)';
    if i==1
        FSegmentsCon = [FSegmentsCon;IDs];
    else
        if ~ismember(FSegmentsCon,[IDs;flip(IDs)],'rows')
            FSegmentsCon = [FSegmentsCon;IDs];
        else
            rm_IDs = [rm_IDs;i];
        end
    end
      
end
Segment(rm_IDs,:) = [];
NSegment = Segment;

% Derive Junction Connectivity
nJunctions = size(FJunctions,1);
FOnEdge = zeros(nJunctions,1);
for i = 1:nJunctions
    r = FJunctions(i,:);
    IsOnEdge1 = find(r==0);
    IsOnEdge2 = find(r==L);
    if ~isempty(IsOnEdge1) || ~isempty(IsOnEdge2)
        FOnEdge(i) = 1;
    end
    AsStart = find(i==FSegmentsCon(:,1));
    AsEnd = find(i==FSegmentsCon(:,2));
    TempSegments=setxor(union(AsStart,AsEnd),intersect(AsStart,AsEnd));
    TempJunctions = unique(FSegmentsCon(TempSegments,:));
    FJunctiosCon{i} =setxor(i,TempJunctions);
    FLocConNum(i) = length(FJunctiosCon{i});  
end

% Inner Junctions connectivity
InnerConNum = FLocConNum(~FOnEdge);