function [X_new] = terminate(SegEnds,X,FaceOut)
%% Code
switch size(FaceOut,2)
    %% segment for which one coordinates out of bound. 
    case 1 
        switch FaceOut
            case 1
                ratio=abs((1-SegEnds(1))/(X(1)-SegEnds(1)));
                X_new=SegEnds+ratio*(X-SegEnds);
            case 2
                ratio=abs((0-SegEnds(1))/(X(1)-SegEnds(1)));
                X_new=SegEnds+ratio*(X-SegEnds);
            case 3
                ratio=abs((1-SegEnds(2))/(X(2)-SegEnds(2)));
                X_new=SegEnds+ratio*(X-SegEnds);
            case 4
                ratio=abs((0-SegEnds(2))/(X(2)-SegEnds(2)));
                X_new=SegEnds+ratio*(X-SegEnds);
            case 5
                ratio=abs((1-SegEnds(3))/(X(3)-SegEnds(3)));
                X_new=SegEnds+ratio*(X-SegEnds);
            case 6
                ratio=abs((0-SegEnds(3))/(X(3)-SegEnds(3)));
                X_new=SegEnds+ratio*(X-SegEnds);         
        end
%% segment for which two coordinates out of bound.       
    case 2 
        %% Face 1 and adjacent edges
        if (FaceOut(1)==1 && FaceOut(2)==3) || (FaceOut(1)==3 && FaceOut(2)==1)
           ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(2))/(X(2)-SegEnds(2));
           ratio=min(ratio1,ratio2);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==1 && FaceOut(2)==4) || (FaceOut(1)==4 && FaceOut(2)==1)
            ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(2))/(X(2)-SegEnds(2));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==1 && FaceOut(2)==5) || (FaceOut(1)==5 && FaceOut(2)==1)
            ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==1 && FaceOut(2)==6) || (FaceOut(1)==6 && FaceOut(2)==1)
            ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        %% Face 2 and adjacent edges
        if (FaceOut(1)==2 && FaceOut(2)==3) || (FaceOut(1)==3 && FaceOut(2)==2)
            ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(2))/(X(2)-SegEnds(2));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==2 && FaceOut(2)==4) || (FaceOut(1)==4 && FaceOut(2)==2)
            ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(2))/(X(2)-SegEnds(2));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==2 && FaceOut(2)==5) || (FaceOut(1)==5 && FaceOut(2)==2)
            ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==2 && FaceOut(2)==6) || (FaceOut(1)==6 && FaceOut(2)==2)
            ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        %% Face 3 and adjacent edges (except those considered before)
        if (FaceOut(1)==3 && FaceOut(2)==5) || (FaceOut(1)==5 && FaceOut(2)==3)
            ratio1=(1-SegEnds(2))/(X(2)-SegEnds(2));ratio2=(1-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==3 && FaceOut(2)==6) || (FaceOut(1)==6 && FaceOut(2)==3)
            ratio1=(1-SegEnds(2))/(X(2)-SegEnds(2));ratio2=(0-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        %% Face 4 and adjacent edges (except those considered before)
        if (FaceOut(1)==4 && FaceOut(2)==5) || (FaceOut(1)==5 && FaceOut(2)==4)
            ratio1=(0-SegEnds(2))/(X(2)-SegEnds(2));ratio2=(1-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
        if (FaceOut(1)==4 && FaceOut(2)==6) || (FaceOut(1)==6 && FaceOut(2)==4)
            ratio1=(0-SegEnds(2))/(X(2)-SegEnds(2));ratio2=(0-SegEnds(3))/(X(3)-SegEnds(3));
            ratio=min(ratio1,ratio2);
            X_new=SegEnds+ratio*(X-SegEnds);
        end
%% NOTE: Above are all possible face pairs. (12 face pairs = 12 edges)       
 %% segment for which three coordinates out of bound.        
    case 3 
        % Face 1,3 and 5
        if sum(ismember([1,3,5],FaceOut))==3
%            disp ('Face 1,3 and 5')
           ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(1-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 1,3 and 6
        if sum(ismember([1,3,6],FaceOut))==3
%            disp ('Face 1,3 and 6')
           ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(0-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 1,4 and 5
        if sum(ismember([1,4,5],FaceOut))==3
%            disp ('Face 1,4 and 5')
           ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(1-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 1,4 and 6
        if sum(ismember([1,4,6],FaceOut))==3
%            disp ('Face 1,4 and 6')
           ratio1=(1-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(0-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 2,3 and 5
        if sum(ismember([2,3,5],FaceOut))==3
%            disp ('Face 2,3 and 5')
           ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(1-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 2,3 and 6
        if sum(ismember([2,3,6],FaceOut))==3
%            disp ('Face 2,3 and 6')
           ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(1-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(0-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 2,4 and 5
        if sum(ismember([2,4,5],FaceOut))==3
%            disp ('Face 2,4 and 5')
           ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(1-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end
        % Face 2,4 and 6
        if sum(ismember([2,4,6],FaceOut))==3
%            disp ('Face 2,4 and 6')
           ratio1=(0-SegEnds(1))/(X(1)-SegEnds(1));ratio2=(0-SegEnds(2))/(X(2)-SegEnds(2));
           ratio3=(0-SegEnds(3))/(X(3)-SegEnds(3));ratio=min(min(ratio1,ratio2),ratio3);
           X_new=SegEnds+ratio*(X-SegEnds);
        end  
end
% disp(X_new)
end

