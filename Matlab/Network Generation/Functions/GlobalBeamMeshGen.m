function [GlobalBeamMesh] = GlobalBeamMeshGen(BeamMesh, dist_tol)
% Receives the meshed beam domains (branches) and generates a global mesh.
% Merges common nodes (aka crosslinks) (JUST the nodes, not the elements).
% For constant stress transfer test (consistency test) merging overlapping
% nodes is not applied.

n_beams = length(BeamMesh);

% Loop beam branches
elem_count = 0;
node_count = 0;
GlobalNodes = [];
GlobalCon = [];
GlobalSParam = [];
GlobalTParam = [];
GlobalSegmentMat = {};


% Different n1, n2 for quadratic beams
if BeamMesh(1).Order==2
    [~,dN,dN2] = ShapeFnc(-1, 2); % Export only for the first node
end

for i = 1:n_beams
    
    % Export Beam branch data
    nB = BeamMesh(i).nB;
    Nodes = BeamMesh(i).Nodes  ;
    Connectivity = BeamMesh(i).Connectivity;  
    
    % Merge matrices
    GlobalNodes = [GlobalNodes;Nodes];
    GlobalCon = [GlobalCon; Connectivity+node_count];
    
    % Update nodecounter
    node_count = node_count+nB;

    % Global Element quantities (are NOT merged, no matter what)
    nel = size(Connectivity,1);
    GlobalSParam = [GlobalSParam;BeamMesh(i).SParam];   
   
end

% Export Data
GlobalBeamMesh.Nodes = GlobalNodes;
GlobalBeamMesh.SParam = GlobalSParam;
GlobalBeamMesh.TParam = GlobalTParam;
GlobalBeamMesh.Connectivity = GlobalCon;
GlobalBeamMesh.SegmentMat=GlobalSegmentMat;



% %% Merge overlapping nodes if requested
% Group2Merge = {};
% counter = 0;
% nBg = size(GlobalNodes,1);
% for i = 1:nBg
%     [~, D] = knnsearch(GlobalNodes(i,:), GlobalNodes);
%     IDX = find(D<dist_tol); %IDX always includes [i]-node
% 
%     if length(IDX)>1 % two or more nodes to be merged
%         counter = counter+1;
%         Group2Merge{counter} = sort(IDX)';
%     end
% end
% 
% % Keep only unique Groups
% a = Group2Merge(:);
% [a1,b,c] = unique(cellfun(@char,a,'un',0));
% lo = histc(c,1:max(c));
% loo = lo(:) > 1;
% UniqueGroups = [a(b(loo)), num2cell(lo(loo))];
% 
% % Replace that in connectivity matrix and Node matrix
% IDs = (1:nBg)';
% rows2erase = [];
% for i = 1:size(UniqueGroups,1)
%     mergedIDs = UniqueGroups{i,1};
%     n_merged = length(mergedIDs);
%     fprintf('Merged %d beam nodes \n',n_merged)
%     for ii = 2:n_merged
%         IDs(mergedIDs(ii)) = mergedIDs(1);
%         rows2erase = [rows2erase,mergedIDs(ii)];
%     end
% end
% % Find Deleted nodes
% IDXs = find(IDs~=(1:nBg)');
% newIDs = IDs;
% newIDs(IDXs) = [];
% 
% for i = 1:length(newIDs)
%     ID2keep = find(newIDs(i)==IDs);
%     IDs(ID2keep)=i;
% end
% 
% %Finally, define new connectivity
% NewGlobalCon = IDs(GlobalCon);
% 
% % Update Global nodes
% NewGlobalNodes = GlobalNodes;
% NewGlobalNodes(IDXs,:) = [];
% 
% % Update Output structure
% GlobalBeamMesh.Nodes = NewGlobalNodes;
% GlobalBeamMesh.Connectivity = NewGlobalCon;
% GlobalBeamMesh.Nodes2Merge = UniqueGroups;
% 
% 
% 
% 
% 
