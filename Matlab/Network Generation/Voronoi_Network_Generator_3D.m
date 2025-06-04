%% Voronoi 3D network Generator
clc
clear
close all
addpath 'Functions/'
addpath 'VonNetworkFunctions/'

% This script requires the Matlab package: matGeom by David Legland

%% INPUT
G = 1.01; % Size of seed generation domain
L = 1; % Characteristic length of network
nseeds = 60 % Number of seeds

% Loop to create multiple networks
rng default
for rng_loop = 1:11

% Use this variable to pick where you save your networks and what the
% variable's name will be. Output is a structure Network that contains the
% fiber nodal coordinates (Network.Segments), along with the parameters
% used to create the network, and other things that are not needed to
% create a mesh and run sims.
network_folder = [num2str(nseeds),'/'];
network_name = ['Case',num2str(rng_loop)];


rm_bool = true; % Boolean to remove Junctions or not
z_target = 3.4;% Coordination number
n_seg_rm = 4; %Random batch removal of segments
CO_length = .01;
round_off_edge_tol = 1e-2;

% Prealigned
lambda = 1;
ref_vec = [0,1,0];

%% GENERATE INITIAL NETWORK
Vcub=[0 0 0; L 0 0; 0 L 0; 0 0 L;
      L L L; 0 L L; L 0 L; L L 0];
seeds=G*rand(nseeds,3)-(G-L)/2;

%% Generate the initial network
Feff = diag([lambda,lambda^(-1/2),lambda^(-1/2)]);
[Ver,Cel,C_tst,ok_cells,seeds_tf]=voronoi3d_cuboid(seeds,Vcub, 12,Feff);

%% Extract the segments
% Initialize
SegmentsRaw = [];
fib_count = 0;
fib_length = [];

figure
hold on
% axis('equal')
% view([-36 27])
% axis off
% scatter3(seeds_tf(:,1),seeds_tf(:,2),seeds_tf(:,3),25, ...
%          'Marker','o','MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','k');
% scatter3(Ver(:,1),Ver(:,2),Ver(:,3),25, ...
%          'Marker','o','MarkerFaceColor',[0 1 0], 'MarkerEdgeColor','r');
%          counter= 0;

for kl = 1:length(ok_cells)
    k = ok_cells(kl);
if ~isempty(Cel{k})
    col=rand(1,3);
    Vk = Ver(Cel{k},:); Fk = convhull(Vk);
    
    [Vk, Fk] = mergeCoplanarFaces(Vk, Fk);
    for i=1:length(Fk)
        patch('Vertices',Vk,'Faces',Fk{i},'FaceColor',col,'FaceAlpha',0.3) 

        
        temp = Fk{i};
                    
        for ii = 1:length(temp)
            r0 = Vk(temp(ii),:);
            if ii<length(temp)
                rf = Vk(temp(ii+1),:);
            else
                rf = Vk(temp(1),:);
            end
            fib = rf-r0;
            OnEdge = false;
            ParallelEdge = find(fib==0);
            if ~isempty(ParallelEdge)
                OnEdge01 = find(r0==1);
                OnEdge00 = find(r0==0);
                OnEdgef1 = find(rf==1);
                OnEdgef0 = find(rf==0);
                if ~isempty(OnEdge01)||~isempty(OnEdge00)||~isempty(OnEdgef1)||~isempty(OnEdgef0)
                    OnEdge = true;
                end

            end
         
            if ~OnEdge
                fib_count = fib_count+1;
                SegmentsRaw =[SegmentsRaw;r0,rf];
            end
        end
    end
    
end
end
xlabel('X');ylabel('Y');zlabel('Z');

%% ORGANIZE NETWORK DATA
% Define Connectivity at Junction
Junctions = unique([SegmentsRaw(:,1:3);SegmentsRaw(:,4:6)],'rows');
SegmentsCon = [];

% Define Unique Segments Connectivity
for i = 1:size(SegmentsRaw,1)
    r0 = SegmentsRaw(i,1:3);
    rf = SegmentsRaw(i,4:6);
    index=ismember(Junctions,[r0;rf],'rows');
    IDs = find(index)';
    if i==1
        SegmentsCon = [SegmentsCon;IDs];
    else
        if ~ismember(SegmentsCon,[IDs;flip(IDs)],'rows')
            SegmentsCon = [SegmentsCon;IDs];
        end
    end
      
end

% Apply Fiber Length Restrictions
n_fibers = size(SegmentsCon,1);
SegmentLength = zeros(n_fibers,1);
for i = 1:n_fibers

    r0 = Junctions(SegmentsCon(i,1),:);
    rf = Junctions(SegmentsCon(i,2),:);

    SegmentLength(i) = norm(rf-r0);
end

% CutOffLength =  find(SegmentLength<mean(SegmentLength)/100);
CutOffLength =  find(SegmentLength < CO_length);
MergSegCon = SegmentsCon(CutOffLength,:);
MergedSegs = [];
JunctionMerged = zeros(size(Junctions,1));
for iel = 1:length(CutOffLength)
    
    HasBeenMerged = find(MergedSegs==iel);

    if isempty(HasBeenMerged)
        i = CutOffLength(iel);
    
        Id1 = SegmentsCon(i,1);
        Id2 = SegmentsCon(i,2);
        [SegsW1,~] = find(MergSegCon==Id1);
        [SegsW2,~] = find(MergSegCon==Id2);
        MergedSegs = [union(MergedSegs,[SegsW1;SegsW2])];
        
        % First repetition of merged
        IDs2Merge = union(unique(MergSegCon(SegsW1,:)),unique(MergSegCon(SegsW2,:)));
        IDs2Merge= IDs2Merge(:);

        %Second rep.
        for ii = 1:length(IDs2Merge)
            [idcs,~] = find(MergSegCon==IDs2Merge(ii));
            tempIDs = unique(MergSegCon(idcs,:));
            PotentialNewIDs = setdiff(tempIDs,IDs2Merge);
            for iii = 1:length(PotentialNewIDs)
                IDs2Merge = [IDs2Merge;PotentialNewIDs(iii)];
                MergedSegs = union(MergedSegs,idcs);
            end
        
        end
                        
        rm = mean(Junctions(IDs2Merge,:));
    
        %Replace the Coords with the average of the merged nodes
        Junctions(IDs2Merge,:) = repmat(rm,[length(IDs2Merge),1]);
        JunctionMerged(IDs2Merge) = 1;
        
        for ii = 1:length(IDs2Merge)
            idcs = (SegmentsCon==IDs2Merge(ii));
            SegmentsCon(idcs) = min(IDs2Merge);
        end
    end
        
end

% Check all Merges have been okay
for iel = 1:length(CutOffLength)
    i = CutOffLength(iel);
    r0 = Junctions(SegmentsCon(i,1),:);
    rf = Junctions(SegmentsCon(i,2),:);
    if norm(r0-rf)~=0
        error('Merge is not working: Cut off length is unreasonably large')
    end
end

SegmentLength(CutOffLength) = 0;

% Derive Junction Connectivity
nJunctions = size(Junctions,1);
OnEdge = zeros(nJunctions,1);
for i = 1:nJunctions
    r = Junctions(i,:);
    IsOnEdge1 = find(r==L);
    IsOnEdge2 = find(r==0);
    if ~isempty(IsOnEdge1) || ~isempty(IsOnEdge2)
        OnEdge(i) = 1;
    end
    AsStart = find(i==SegmentsCon(:,1));
    AsEnd = find(i==SegmentsCon(:,2));
    TempSegments=setxor(union(AsStart,AsEnd),intersect(AsStart,AsEnd));
    TempJunctions = unique(SegmentsCon(TempSegments,:));
    if length(TempJunctions)>2
        
    end
    JunctiosCon{i} =setxor(i,TempJunctions);
    LocConNum(i) = length(JunctiosCon{i});
   
end

% Write segment with length control (skip collapsed elements)
Segment = [];
for i=1:size(SegmentsCon,1)
    
    r0 = Junctions(SegmentsCon(i,1),:);
    rf = Junctions(SegmentsCon(i,2),:);
    
    if SegmentsCon(i,1)~=SegmentsCon(i,2)
        Edge2Edge = OnEdge(SegmentsCon(i,1)) && OnEdge(SegmentsCon(i,2));
        if ~Edge2Edge
            Segment = [Segment;r0,rf];
        end
  
    end
end


%% DEFINE NEW Network Matrices
% Apply Round off corrections
Segment(Segment<round_off_edge_tol)=0;
Segment(Segment>L-round_off_edge_tol)=L;

%Define new network matrices
[Segment,FJunctions,FSegmentsCon,FJunctiosCon,FLocConNum,InnerConNum,FFibLength,theta,FOnEdge] = arrange_network(Segment,L,ref_vec);
plot_network(Segment,'b',2)

%% REMOVE SEGMENTS FOR TARGET CONNECTIVITY MATRIX

% % If a segment connects twon Junctions with connectivity>4 remove it
% rm_seg =[];
% for i = 1:size(FSegmentsCon,1)
%     ID1 = FSegmentsCon(i,1) ;
%     ID2 = FSegmentsCon(i,2) ;
%     if FLocConNum(ID1)>4 || FLocConNum(ID2)>4
%         rm_seg =[rm_seg;i];
%     end
% end
% plot_network(Segment(rm_seg,:),'g',2)
% Segment(rm_seg,:)=[];
% [FJunctions,FSegmentsCon,FJunctiosCon,FLocConNum,InnerConNum,FFibLength,theta,FOnEdge] = arrange_network(Segment,L,ref_vec);

if rm_bool
    z_temp = mean(InnerConNum);
    
    stop_bool = false;
    
    iter = 0;
    while ~stop_bool 
        iter = iter + 1;
        fprintf('Segment removal iteration ...#%d \n', iter)
        if z_temp<=z_target
            break
        end
        % Find candidate segments for removal
        rm_seg =[];
        for i = 1:size(FSegmentsCon,1)
            ID1 = FSegmentsCon(i,1) ;
            ID2 = FSegmentsCon(i,2) ;
            if FLocConNum(ID1)>=4 && FLocConNum(ID2)>=4
                rm_seg =[rm_seg;i];
            end
        end
    
        if isempty(rm_seg) || length(rm_seg)<n_seg_rm
            break
        end
        
        ID_rm = randi(length(rm_seg),n_seg_rm,1);
%         plot_network(Segment(rm_seg(ID_rm),:),'r',2)
        Segment(rm_seg(ID_rm),:) = [];
        
    
        [Segment,FJunctions,FSegmentsCon,FJunctiosCon,FLocConNum,InnerConNum,FFibLength,theta,FOnEdge] = arrange_network(Segment,L,ref_vec);
        z_temp = mean(InnerConNum) %#ok<*NOPTS>
    end

end
%% REMOVE JUNCTIONS: DISCONNECTED COMPONENTS AND DANGLING ENDS
Junctions2Delet = ~FOnEdge & (FLocConNum<=2)';

count1 = 0 ;
while sum(Junctions2Delet)~=0

%     %--------------------plot---------------------------------------------
%     count1 = count1+1;
%     plot_network(Segment,'b',1)
%     hold on
%     plot3(FJunctions(Junctions2Delet,1),FJunctions(Junctions2Delet,2),FJunctions(Junctions2Delet,3),'r.','markersize',25)
%     %----------------------------------------------------------------------
    
    % Delete or merge Junction    
    tempIDs = find(Junctions2Delet);
    for i=1:length(tempIDs)
        tempCon = FJunctiosCon{tempIDs(i)};
        if length(tempCon)>1
           if ~ismember(FSegmentsCon,[tempCon';flip(tempCon')],'rows')
               %the segment does not exist already.
               % Instead of deleting the junction, connect directly the two
               % otherwise inderectly connected junctions, by merging:
               % Also, check length requirement
               newSegmL = norm(FJunctions(tempCon(1),:)-FJunctions(tempCon(2),:));
               if newSegmL>CO_length
                    FSegmentsCon(FSegmentsCon==tempIDs(i))=tempCon(1);
               else
                    FJunctions(tempIDs(i),:) = nan;
               end
           else
               % Segment exists already, so just delete the dangling end
               FJunctions(tempIDs(i),:) = nan;
           end
        else
            %Just one hanging fiber
            FJunctions(tempIDs(i),:) = nan;
        end
    end

    % Write new segment matrix
    Segment = [];
    for i=1:size(FSegmentsCon,1)
        
        r0 = FJunctions(FSegmentsCon(i,1),:);
        rf = FJunctions(FSegmentsCon(i,2),:);
        
        if ~isnan(r0(1)) && ~isnan(rf(1)) && FSegmentsCon(i,1)~=FSegmentsCon(i,2)
            Edge2Edge = FOnEdge(FSegmentsCon(i,1)) && FOnEdge(FSegmentsCon(i,2));
            if ~Edge2Edge
                Segment = [Segment;r0,rf];
            end
        end
    end
    
    [Segment,FJunctions,FSegmentsCon,FJunctiosCon,FLocConNum,~,FFibLength,~,FOnEdge] = arrange_network(Segment,L,ref_vec);
    Junctions2Delet = ~FOnEdge & (FLocConNum<=2)';
end

%% Export one final time
[Segment,Junctions,SegmentsCon,JunctiosCon,LocConNum,InnerConNum,FibLength,theta,OnEdge] = arrange_network(Segment,L,ref_vec);
% CutOffLength =  find(FibLength < 0.1);
% Segment(CutOffLength,:) = [];
Network.Segment = Segment;
Network.Junction = Junctions;
Network.JunctiosCon = JunctiosCon;
Network.LocConNum = LocConNum;
Network.InnerConNum = InnerConNum;
Network.FibLength = FibLength;
Network.theta = theta;
Network.OnEdge = OnEdge;
% Options
Network.NSeeds = nseeds;
Network.L = L;
% Network.G = G;
Network.rm_bool = rm_bool; % Boolean to remove Junctions or not
Network.z_target = z_target;
Network.CO_length = CO_length;
Network.lambda =lambda;
Network.ref_vec = ref_vec;
Network.round_off_edge_tol = round_off_edge_tol;

% Adjust Segment and Junction prior to saving Network variable
% Remove small fibers and strays - note, fibers smaller than an element
% length will cause problems with convergence
Network.Segment = remove_small_fibers(Network.Segment);
Network.Segment = remove_stray_fibers(Network.Segment);
% Adjust Junction after removing small and stray fibers
Network = adjust_network(Network);

Segment = Network.Segment; Junction = Network.Junction;

% Save in both folders to have in case undulation is desired
if ~isfolder(['../Networks/Data/3D/Straight/',network_folder])
    mkdir(['../Networks/Data/3D/Straight/',network_folder]); end 

if ~isfolder(['../Networks/Data/3D/Undulated/',network_folder])
    mkdir(['../Networks/Data/3D/Undulated/',network_folder]); end 

save(['../Networks/Data/3D/Straight/',network_folder,network_name,'.mat'],'Segment','Junction')
save(['../Networks/Data/3D/Undulated/',network_folder,network_name,'.mat'],'Segment','Junction')

% Final Statistics
zfinal = mean(InnerConNum);
lc = mean(FFibLength);

fprintf('Number of fibers is %d. Network mean fiber length is %f. Final connectivity is z = %f\n',length(Segment),lc,zfinal)

end


