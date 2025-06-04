%% CUBE COORDINATE MAPPING TRIALS
function C = BCC_lat(Nx,Ny,Nz,Lx,Ly,Lz)
O = [0 0 0];% origin

%% COORDINATES FOR BODY-CENTERED CUBIC (BCC) ELEMENT

nodes = [O; O(1)+Lx, O(2), O(3); O(1)+Lx, O(2)+Ly, O(3); O(1), O(2)+Ly, O(3);...
    O(1), O(2), O(3)+Lz; O(1)+Lx, O(2), O(3)+Lz; O(1)+Lx, O(2)+Ly, O(3)+Lz;...
    O(1), O(2)+Ly, O(3)+Lz; O(1)+Lx/2, O(2)+Ly/2, O(3)+Lz/2];

no_nodes = length(nodes); 
%% LATTICE NODES

% cell replication - produce every node in lattice (including duplicates)

for nix=1:Nx-1  %Modified to include variable number of nodes for each cell type
    for i=1:no_nodes
            nodes(no_nodes+no_nodes*(nix-1)+i,:)=[nodes(i,1)+nix*Lx,nodes(i,2),nodes(i,3)];
    end
end
for niy=1:Ny-1
    for i=1:no_nodes*Nx
            nodes(no_nodes*Nx+no_nodes*Nx*(niy-1)+i,:)=[nodes(i,1),nodes(i,2)+niy*Ly,nodes(i,3)];
    end
end
for niz=1:Nz-1
    for i=1:no_nodes*Nx*Ny
            nodes(no_nodes*Nx*Ny+no_nodes*Nx*Ny*(niz-1)+i,:)=[nodes(i,1),nodes(i,2),nodes(i,3)+niz*Lz];
    end
end
%% Remove duplicate nodes
[C,~,~]=unique(nodes,'rows');
% %% SCATTER NODES
% 
% figure
% scatter3(C(:,1),C(:,2),C(:,3),'rx');
% hold on
% C1 = C(:,1); C2 = C(:,2); C3 = C(:,3);
% 
% for k=1:length(C)
%     text(C1(k),C2(k),C3(k),num2str(k))
% end
% 
% daspect([1 1 1])
% xlabel('xaxis'); ylabel('yaxis'); zlabel('zaxis');
% % xlim([-0.5*L,Nx])
% % ylim([-0.5*L,Ny])
% % zlim([-0.5*L,2])