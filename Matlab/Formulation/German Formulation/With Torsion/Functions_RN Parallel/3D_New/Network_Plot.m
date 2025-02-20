function [] = Network_Plot( P,Q,Patch,Patch_Bs)
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\C_files');
addpath('C:\Users\smm5969\Box Sync\Fibrin Modeling\German Formulation\nurbs_toolbox');
% %% Plot undeformed shape
figure
% pts=20;
% for i=1:size(Patch,2)
% % Find NURBS shape fun.
% order=Patch(i).order;kv=Patch(i).kv;w=Patch(i).wts;n=Patch(i).num_cpt;
% nos=size(n,2);
% for k=1:nos
%     for j=1:pts+1
%         xi=(1/pts)*(j-1);
%         [N(j,k),dN(j,k)]=NURBSbasis(k,order,xi,kv,w);
%     end
% end
% P=Patch(i).cpt;
% P=P(1:3,:);
% L0=zeros(pts+1,3);
% for j=1:nos
%     L0(:,1)=N(:,j)*P(1,j)+L0(:,1);L0(:,2)=N(:,j)*P(2,j)+L0(:,2);L0(:,3)=N(:,j)*P(3,j)+L0(:,3);
% end
% if i<=size(Patch,2)-size(Patch_Bs,2)
%     plot(L0(:,1),L0(:,2),'b-','LineWidth',4);hold on
% else
%     plot(L0(:,1),L0(:,2),'g-','LineWidth',4);hold on
% end
% end
% Plot deformed shape
pts=20;
for i=1:size(Patch,2)
% Find NURBS shape fun.
order=Patch(i).order;kv=Patch(i).kv;w=Patch(i).wts;n=Patch(i).num_cpt;
nos=size(n,2);
for k=1:nos
    for j=1:pts+1
        xi=(1/pts)*(j-1);
        [N(j,k),dN(j,k)]=NURBSbasis(k,order,xi,kv,w);
    end
end
P=Q(:,n);
P=P(1:3,:);
L0=zeros(pts+1,3);
for j=1:nos
    L0(:,1)=N(:,j)*P(1,j)+L0(:,1);L0(:,2)=N(:,j)*P(2,j)+L0(:,2);L0(:,3)=N(:,j)*P(3,j)+L0(:,3);
end
if i<=size(Patch,2)-size(Patch_Bs,2)
    plot(L0(:,1),L0(:,2),'r','LineWidth',4);hold on
else
    plot(L0(:,1),L0(:,2),'r','LineWidth',4);hold on
end
end
% title('Undeformed (Black) and Deformed Shape(Blue)')
xlabel('X');ylabel('Y');
end









