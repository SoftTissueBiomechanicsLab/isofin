function [ R ] = Network_R_New(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,DOF_Bs,D_DOF,F_DOF)
%% Code
%% Modify Q as per Dirichlet's BC
% % % for i=1:size(P,2)
% % %     for j=1:4
% % %         if D_DOF(i,j)~=Inf
% % %         % Modify Q
% % %            Q(j,i)=P(j,i)+D_DOF(i,j);
% % %         end
% % %     end
% % % end
size_Patch=size(Patch);size_Patch_Bs=size(Patch_Bs);
%% Form K and R without bending strip for the whole network.
% Store tangent stiffness and residue for all patches.
K_cell={};R_cell={}; % Cells to store all K_patch and R_patch.
parfor i=1:size_Patch(2) %% PARFOR LOOP
    % Find K_patch and R_patch.
    n=Patch(i).num_cpt;order=Patch(i).order;kv=Patch(i).kv;wts=Patch(i).wts;
    nos=size(n,2);Nele=nos-order;p=P(:,n);q=Q(:,n);
    End_pt=zeros(Nele,2);
    for j=1:Nele
        End_pt(j,1:2)=[(j-1)/Nele,j/Nele];
    end
    Ele2Cp=zeros(Nele,order+1);
    for j=1:Nele
        for k=1:order+1
            Ele2Cp(j,k)=j+k-1;
        end
    end
    key_K=1;dof=zeros(nos,4);
    for j=1:nos
        dof(j,:)=[key_K,key_K+1,key_K+2,key_K+3];key_K=key_K+4;
    end
    a0=transpose(A0(:,3*(i-1)+1:3*(i-1)+3));ngp=order+1;
    if i<=(size_Patch(2)-size_Patch_Bs(2))% Distinguishing Between Arcs and fibers.
%         disp('Stiffness Contribution from fiber')
         [ R_patch ] = PatchN_R_New(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat,dof,a0,ngp);
    else
%         disp('Stiffness Contribution from arc')
        [ R_patch ] = PatchN_R_New(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Arcs,dof,a0,ngp);
    end
    
    %% Add K-patch and R_patch to the respective cells
% % %     K_cell{1,i}=K_patch;
    R_cell{1,i}=R_patch;
end
% Assemble K and R from the respective cells.
K=zeros(4*size(P,2));R=zeros(4*size(P,2),1);
for i=1:size_Patch(2) %% FOR LOOP
   n=Patch(i).num_cpt;nos=size(n,2);
   key_K=1;
   for j=1:nos
        dof(j,:)=[key_K,key_K+1,key_K+2,key_K+3];key_K=key_K+4;
   end
% Get K_patch and R_patch for ith patch.
% % %    K_patch=K_cell{1,i};
   R_patch=R_cell{1,i};
% % %    K_patch=full(K_patch);
   R_patch=full(R_patch);
   % Add to K and R
    for r=1:nos
        for r1=1:4
            a=n(1,r);
            R(DOF(a,r1),1)=R(DOF(a,r1),1)+R_patch(dof(r,r1),1);
% % %             for s=1:nos
% % %                 for s1=1:4
% % %                     b=n(1,s);
% % %                     if dof(r,r1)<=dof(s,s1)
% % %                        K(DOF(a,r1),DOF(b,s1))=K(DOF(a,r1),DOF(b,s1))+K_patch(dof(r,r1),dof(s,s1));
% % %                        if DOF(a,r1)~=DOF(b,s1)
% % %                           K(DOF(b,s1),DOF(a,r1))=K(DOF(b,s1),DOF(a,r1))+K_patch(dof(r,r1),dof(s,s1));
% % %                        end
% % %                     end
% % %                 end
% % %             end
        end
    end
end
% disp('Det of K after assembling patches')
% det(K)


%% Add bending strip
% Store tangent stiffness and residue for all bending strips.
K_cell={};R_cell={}; % Cells to store all K_bs and R_bs.
parfor i=1:size_Patch_Bs(2) %% PARFOR LOOP
    % Find K_bs and R_bs.
    n=Patch_Bs(i).num_cpt;order=Patch_Bs(i).order;kv=Patch_Bs(i).kv;wts=Patch_Bs(i).wts;
    nos=size(n,2);Nele=nos-order;p=P(:,n);q=Q(:,n);
    End_pt=zeros(Nele,2);
    for j=1:Nele
        End_pt(j,1:2)=[(j-1)/Nele,j/Nele];
    end
    Ele2Cp=zeros(Nele,order+1);
    for j=1:Nele
        for k=1:order+1
            Ele2Cp(j,k)=j+(k-1);
        end
    end
    key_K=1;dof=zeros(nos,4);
    for j=1:nos
        dof(j,:)=[key_K,key_K+1,key_K+2,key_K+3];key_K=key_K+4;
    end
    a0=transpose(A0_Bs(:,3*(i-1)+1:3*(i-1)+3));ngp=order+1;
    [ R_bs ] = PatchN_R_bs_New(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Bs,dof,a0,ngp);
    %% Add K_bs and R_bs to the respective cells
%     K_cell{1,i}=K_bs;
    R_cell{1,i}=R_bs;
end
% Add bending contribution to  K and R from the respective cells.
for i=1:size_Patch_Bs(2) % FOR LOOP
   n=Patch_Bs(i).num_cpt;nos=size(n,2);
   key_K=1;
   for j=1:nos
        dof(j,:)=[key_K,key_K+1,key_K+2,key_K+3];key_K=key_K+4;
   end
% Get K_patch and R_patch for ith patch.
% % %    K_bs=K_cell{1,i};
   R_bs=R_cell{1,i};
% % %    K_bs=full(K_bs);
   R_bs=full(R_bs);
   % Add to K and R
    for r=1:nos
        for r1=1:4
            a=n(1,r);
            R(DOF(a,r1),1)=R(DOF(a,r1),1)+R_bs(dof(r,r1),1);
% % %             for s=1:nos
% % %                 for s1=1:4
% % %                     b=n(1,s);
% % %                     if dof(r,r1)<=dof(s,s1)
% % %                        K(DOF(a,r1),DOF(b,s1))=K(DOF(a,r1),DOF(b,s1))+K_bs(dof(r,r1),dof(s,s1));
% % %                        if DOF(a,r1)~=DOF(b,s1)
% % %                           K(DOF(b,s1),DOF(a,r1))=K(DOF(b,s1),DOF(a,r1))+K_bs(dof(r,r1),dof(s,s1));
% % %                        end
% % %                     end
% % %                 end
% % %             end
        end
    end
end
% disp('Det of K after assembling bending strips')
% det(K)

%% Apply Boundary Conditions
% % Apply Neumann BC (Force BC)
% for i=1:size(P,2)
%     for j=1:4
%         % Modify R
%         key_K=DOF(i,j);
%         R(key_K,1)= R(key_K,1)-F_DOF(i,j);
%     end
%     
% end
% % disp('Det of K after Neumann BC')
% % det(K)
% % Apply Dirichlet's BC (Displacement BC)
% R1=R;
% for i=1:size(P,2)
%     for j=1:4
%         if D_DOF(i,j)~=Inf
%             key_K=DOF(i,j);
%         % Modify K
% % % %            K(key_K,:)=0; K(key_K,key_K)=1;
% %            K(key,key)=1;
%         % Modify R
%            R(key_K,:)=0;
%         end
%     end
% end
% disp('Det of K after Dirichlet BC')
% det(K)
end






