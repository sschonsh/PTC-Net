function [Surf] = DefineUpsampleInterp(Surf)

for level = 1:length(Surf.Sample)-1
    InterpUp = sparse(Surf.nPts,Surf.nPts);
    for iter = 1:length(Surf.Sample{level})
        pt = Surf.Sample{level}(iter);
        if ismember(pt,Surf.Sample{level+1})
            %keep value
            InterpUp(pt,pt) = 1;
        else
            %find pts
            numLocals = 10;
            NearestUpsampled = knnsearch(Surf.pts(Surf.Sample{level+1},:),Surf.pts(pt,:),'K',numLocals );
            %project onto tangent around pt
            LocalCords = Surf.pts(NearestUpsampled,:)-repmat(Surf.pts(pt,:),numLocals ,1);
            L1 = sum(LocalCords.*repmat(Surf.VF{1}(pt,:),numLocals,1),2);
            L2 = sum(LocalCords.*repmat(cross(Surf.VF{1}(pt,:),Surf.normals(pt,:)),numLocals,1),2);
            
            %find triangulation
            TR = delaunayTriangulation(L1,L2);
            for TRGidx = 1:length(TR.ConnectivityList)
                Coords = cartesianToBarycentric(TR,TRGidx,[0,0]);
                if min(Coords) >= 0
                    %correct weights
                   InterpUp(pt,Surf.Sample{level+1}(NearestUpsampled(TR(TRGidx,:)))) = Coords;%THIRD ARGUEMENT 
                   break
                end   
            end

        end      
    end
    Surf.InterpUp{level} = InterpUp;
end