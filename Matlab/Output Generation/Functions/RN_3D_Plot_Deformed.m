function [] = RN_3D_Plot_Deformed(Q,Patch,Patch_Bs)
%% Code
addpath('../Formulation/German Formulation/C_files')
addpath('../Formulation/German Formulation/nurbs_toolbox')
figure
%% Plot deformed shape
pts=20;
for i=1:size(Patch,2)
    % Find NURBS shape fun.
    order=Patch(i).order;kv=Patch(i).kv;w=Patch(i).wts;n=Patch(i).num_cpt;
    nos=size(n,2);
    N=[];dN=[];
    for k=1:nos
        for j=1:pts+1
            xi=(1/pts)*(j-1);
            [N(j,k),dN(j,k)]=NURBSbasis(k,order,xi,kv,w);
        end
    end
    q=Q(1:3,n);
    L0=zeros(pts+1,3);
    for j=1:nos
        L0(:,1)=N(:,j)*q(1,j)+L0(:,1);L0(:,2)=N(:,j)*q(2,j)+L0(:,2);L0(:,3)=N(:,j)*q(3,j)+L0(:,3);
    end
    % if i<=size(Patch,2)-size(Patch_Bs,2)
    %     plot3(L0(:,1),L0(:,2),L0(:,3),'g','LineWidth',3);hold on
    % else
    %     plot3(L0(:,1),L0(:,2),L0(:,3),'g','LineWidth',3);hold on
    % end
    plot3(L0(:,1),L0(:,2),L0(:,3),'g','LineWidth',3);hold on
end
% title('Undeformed (Black) and Deformed Shape(Blue)')
xlabel('X');ylabel('Y');
end

