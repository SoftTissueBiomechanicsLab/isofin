function [ Reactions ] = Network_Residue(P,Q,Patch,Patch_Bs,A0,A0_Bs,Mat,Mat_Bs,DOF,DOF_Bs,D_DOF,F_DOF)
%% Form R without bending strip
R=zeros(4*size(P,2),1);
parfor i=1:size(Patch,2)
    % Find R_patch.
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
    key=1;dof=zeros(nos,4);
    for j=1:nos
        dof(j,:)=[key,key+1,key+2,key+3];key=key+4;
    end
    a0=transpose(A0(:,3*(i-1)+1:3*(i-1)+3));ngp=order+1;
    if i<=(size(Patch,2)-size(Patch_Bs,2))
         [ R_patch ] = PatchN_R(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat,dof,a0,ngp);
    else
        [ R_patch ] = PatchN_R(Nele,p,q,End_pt,Ele2Cp,kv,order,wts,Mat_Bs,dof,a0,ngp);
    end
    % Add to R
    R_temp=zeros(4*size(P,2),1);
    for r=1:nos
        for r1=1:4
            a=n(1,r);
            R_temp(DOF(a,r1),1)=R_temp(DOF(a,r1),1)+R_patch(dof(r,r1),1);
%             for s=1:nos
%                 for s1=1:4
%                     b=n(1,s);
%                     K(DOF(a,r1),DOF(b,s1))=K(DOF(a,r1),DOF(b,s1))+K_patch(dof(r,r1),dof(s,s1));
%                 end
%             end
        end
    end
    R=R+R_temp;
end

%% Calculate Reactions (fiber ends)
Reactions=zeros(size(D_DOF));
for i=1:size(D_DOF,1)
    for j=1:4
        if D_DOF(i,j)==0
            Reactions(i,j)=R(DOF(i,j),1);
        end
    end
end
end






