function [Surf] = DefineLocalInfoWithRestriction(Surf)
% %
% path(path, 'geodesics/toolbox_fast_marching/');
% path(path, 'toolbox_fast_marching/toolbox/');
% path(path, 'geodesics/toolbox_graph/');
% path(path, 'geodesics/toolbox_graph/off/');
% path(path, 'geodesics/meshes/');

%Triangluation of Temlate
TR = delaunayTriangulation(Surf.Template(:,1),Surf.Template(:,2));

for j = 1:length(Surf.pts)
    %Find Geodesics
    [~, EucDis] = rangesearch(Surf.pts,Surf.pts(j,:),Surf.TempRad*1.25);%find local Euc ball
    [~, pPtsIdx] = max(cell2mat(EucDis));%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    options.end_points = pPtsIdx;%set end point as max in local ball
    center = j;%Surf.pts(j,:);
    [D,~,~]  = perform_fast_marching_mesh(Surf.pts, Surf.trg, center,options);%fast marching
    Surf.interps{j} = find(D < Surf.TempRad*.95);%only keep ones with small local ball
%     Surf.interps{j} = cell2mat(possiblePts)';
    numLocals = length(Surf.interps{j});
    
    %Define Sparse restriction matrix
    R = sparse(Surf.nPts,length(Surf.interps{j}));
    for spPoints = 1:length(Surf.interps{j})
        R(Surf.interps{j}(spPoints),spPoints) = 1;
    end
    Surf.R{j} = R;

    %Compute LCoeffs on each VF
    for k = 1:Surf.nVFs
        %define local coordinates for inter ps
        LocalCords = Surf.pts(Surf.interps{j},:)-repmat(Surf.pts(j,:),numLocals,1);
        Surf.LCoeff{k,j,1} = sum(LocalCords.*repmat(Surf.VF{k}(j,:),numLocals,1),2);
        Surf.LCoeff{k,j,2} = sum(LocalCords.*...
            repmat(cross(Surf.VF{k}(j,:),Surf.normals(j,:)),numLocals,1),2);
        
        %Compute Weights
        W{k,j} = sparse(length(Surf.interps{j}),length(Surf.Template));
        C = squeeze(cell2mat(Surf.LCoeff(k,j,:)));
        for PTidx = 1:length(Surf.LCoeff{k,j,1})
            TRGidx = 1;
            Positive = 0;
            PT = C(PTidx,:);
            %Loop through triangles untill coords are positive
            for TRGidx = 1:length(TR.ConnectivityList)
                Coords = cartesianToBarycentric(TR,TRGidx,PT);
                if min(Coords) >= 0
                    W{k,j}((PTidx),TR(TRGidx,:)) = Coords;%THIRD ARGUEMENT 
                    break
                end   
            end
        end
    end
end
Surf.W = W;
