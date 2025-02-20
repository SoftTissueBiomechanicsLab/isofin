function [K,R] = ParallelAssemble_KR(K,R,DOF,K_patch,R_patch,dof,nos,n)
%% Code
% Add to K and R
for r=1:nos
    for r1=1:4
        a=n(1,r);
        R(DOF(a,r1),1)=R(DOF(a,r1),1)+R_patch(dof(r,r1),1);
        for s=1:nos
            for s1=1:4
                b=n(1,s);
                K(DOF(a,r1),DOF(b,s1))=K(DOF(a,r1),DOF(b,s1))+K_patch(dof(r,r1),dof(s,s1));
            end
        end
    end
end
end

