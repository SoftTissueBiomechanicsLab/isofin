function [ R ] = Network_R(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Arcs,Mat_Bs,DOF,DOF_Bs,D_DOF,F_DOF)
%% Code
%% Modify Q as per Dirichlet's BC
for i=1:size(P,2)
    for j=1:4
        if D_DOF(i,j)~=Inf
        % Modify Q
           Q(j,i)=P(j,i)+D_DOF(i,j);
        end
    end
end
%% Form K and R without bending strip
% % K=zeros(4*size(P,2));
R=zeros(4*size(P,2),1);
for i=1:size(Patch,2)
    % Find K_patch and R_patch.
    n=Patch(i).num_cpt;order=Patch(i).order;kv=Patch(i).kv;wts=Patch(i).wts;
    nos=size(n,2);Nele=nos-order;p=P(:,n);q=Q(:,n);
    for j=1:Nele
        End_pt(j,1:2)=[(j-1)/Nele,j/Nele];
    end
    for j=1:Nele
        for k=1:order+1
            Ele2Cp(j,k)=j+k-1;
        end
    end
    key=1;
    for j=1:nos
        dof(j,:)=[key,key+1,key+2,key+3];key=key+4;
    end
    a0=transpose(A0(:,3*(i-1)+1:3*(i-1)+3));ngp=order+1;
    if i<=(size(Patch,2)-size(Patch_Bs,2))% Distinguishing Between Arcs and fibers.
         [ R_patch ] = PatchN_R(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat,dof,a0,ngp);
    else
        [ R_patch ] = PatchN_R(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Arcs,dof,a0,ngp);
    end
    % Add to K and R
    for r=1:nos
        for r1=1:4
            a=n(1,r);
            R(DOF(a,r1),1)=R(DOF(a,r1),1)+R_patch(dof(r,r1),1);
% %             for s=1:nos
% %                 for s1=1:4
% %                     b=n(1,s);
% %                     K(DOF(a,r1),DOF(b,s1))=K(DOF(a,r1),DOF(b,s1))+K_patch(dof(r,r1),dof(s,s1));
% %                 end
% %             end
        end
    end
end
%% Add bending strip
for i=1:size(Patch_Bs,2)
    % Find K_bs and R_bs.
    n=Patch_Bs(i).num_cpt;order=Patch_Bs(i).order;kv=Patch_Bs(i).kv;wts=Patch_Bs(i).wts;
    nos=size(n,2);Nele=nos-order;p=P(:,n);q=Q(:,n);
    for j=1:Nele
        End_pt(j,1:2)=[(j-1)/Nele,j/Nele];
    end
    for j=1:Nele
        for k=1:order+1
            Ele2Cp(j,k)=j+(k-1);
        end
    end
    key=1;
    for j=1:nos
        dof(j,:)=[key,key+1,key+2,key+3];key=key+4;
    end
    a0=transpose(A0_Bs(:,3*(i-1)+1:3*(i-1)+3));ngp=order+1;
    [ R_bs ] = PatchN_R_bs(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Bs,dof,a0,ngp);
    % Add to K and R
    for r=1:nos
        for r1=1:4
            a=n(1,r);
            R(DOF_Bs(a,r1),1)=R(DOF_Bs(a,r1),1)+R_bs(dof(r,r1),1);
% %             for s=1:nos
% %                 for s1=1:4
% %                     b=n(1,s);
% %                     K(DOF_Bs(a,r1),DOF_Bs(b,s1))=K(DOF_Bs(a,r1),DOF_Bs(b,s1))+K_bs(dof(r,r1),dof(s,s1));
% %                 end
% %             end
        end
    end
end
%% Apply Boundary Conditions
%% Apply Neumann BC (Force BC)
for i=1:size(P,2)
    for j=1:4
        % Modify R
        key=DOF(i,j);
        R(key,1)= R(key,1)-F_DOF(i,j);
    end
end
%% Apply Dirichlet's BC (Displacement BC)
for i=1:size(P,2)
    for j=1:4
        if D_DOF(i,j)~=Inf
            key=DOF(i,j);
        % Modify K
%            K(key,:)=0;K(key,key)=1;
        % Modify R
           R(key,:)=0;
        end
    end
end

end






