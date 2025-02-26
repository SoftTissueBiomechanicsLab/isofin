function [X,EPS,KAPPA,ME,BE,TE,SE] = RN_3D_ParaView(P,Q,Patch,A0,Mat,ele_size)
%% Code
addpath('../Formulation/German Formulation/C_files')
addpath('../Formulation/German Formulation/nurbs_toolbox')

X=cell(1,size(Patch,2));
EPS=cell(1,size(Patch,2));
KAPPA=cell(1,size(Patch,2));
ME=cell(1,size(Patch,2));
BE=cell(1,size(Patch,2));
TE=cell(1,size(Patch,2));
SE=cell(1,size(Patch,2));

for i=1:size(Patch,2)
% Get patch information.
    order=Patch(i).order;kv=Patch(i).kv;
    w=Patch(i).wts;
% n gives global number of the control points that belong to ith patch.
    n=Patch(i).num_cpt;
    nos=size(n,2);ngp=order+1;
% Calculate number of elements on the patch based on the given element size.
    nele=ceil((1/ele_size)*norm(P(1:3,n(1))-P(1:3,n(nos))));
    npts=nele+1;
% Find NURBS shape fun.
    N=zeros(npts+1,nos); dN=N;
    for k=1:nos
        for j=1:npts+1
            xi=(1/npts)*(j-1);
            % Call NURBSbasis.mexw64
            [N(j,k),dN(j,k)]=NURBSbasis(k,order,xi,kv,w); 
        end
    end
    %% Calculate location(coordinates) of the patch.
    p=P(:,n);q=Q(:,n); 
    L0=zeros(npts+1,4);
    for j=1:nos
        L0(:,1)=N(:,j)*q(1,j)+L0(:,1);L0(:,2)=N(:,j)*q(2,j)+L0(:,2);L0(:,3)=N(:,j)*q(3,j)+L0(:,3);L0(:,4)=N(:,j)*q(4,j)+L0(:,4);
    end
    X{i}=L0;
    %% Calculate strain,curvature and energies
    End_pt=zeros(nele,2);
    for j=1:nele
        End_pt(j,1:2)=[(j-1)/nele,j/nele];
    end
    Ele2Cp=zeros(nele,order+1);
    for j=1:nele
        for k=1:order+1
            Ele2Cp(j,k)=j+k-1;
        end
    end
    a0=A0(:,3*(i-1)+1:3*(i-1)+3);
    
    % Calculate and store patch energies
    [ Eps_Patch,Kappa_Patch,ME_Patch,BE_Patch,TE_Patch,SE_Patch ] = Patch_Energy(nele,p,q,End_pt,Ele2Cp,kv,order,w,Mat,a0,ngp);
    
    EPS{i}=Eps_Patch;
    KAPPA{i}=Kappa_Patch;
    ME{i}=ME_Patch;
    BE{i}=BE_Patch;
    TE{i}=TE_Patch;
    SE{i}=SE_Patch;
end

