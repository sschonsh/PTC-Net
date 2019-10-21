function Fun = SpotTranport(K,Surf,VFidx,location,scale,rot)

%Find Correct interps
[possiblePts, ~] = rangesearch(Surf.pts,Surf.pts(location,:),...
    scale*Surf.TempRad*1.1);%find local Euc ball
NewInterps = cell2mat(possiblePts)';
numLocals = length(NewInterps);

%Calculate rotation
CofB = [Surf.VF{VFidx}(location,:);...
    cross(Surf.VF{VFidx}(location,:),Surf.normals(location,:));...
    Surf.normals(location,:)];
Ang  = rot;
Rot  = [cos(Ang), -sin(Ang), 0;...
        sin(Ang),  cos(Ang), 0;...
        0       , 0        , 0;];
B1 = CofB'*Rot*[1,0,0]'; 
B2 = CofB'*Rot*[0,1,0]';



%Find Local Coefficents
LocalCords = Surf.pts(NewInterps,:)-repmat(Surf.pts(location,:),numLocals,1);
LCoeff{1} = sum(LocalCords.*repmat(B1',numLocals,1),2);
LCoeff{2} = sum(LocalCords.*repmat(B2',numLocals,1),2);

%Triangluation
TR = delaunayTriangulation(scale*Surf.Template(:,1),scale*Surf.Template(:,2));        

%Compute Weights
W = sparse(Surf.nPts,length(Surf.Template));
C = squeeze(cell2mat(LCoeff(:)'));
for PTidx = 1:length(LCoeff{1})
    TRGidx = 1;
    Positive = 0;
    PT = C(PTidx,:);
    %Loop through triangles untill coords are positive
    for TRGidx = 1:length(TR.ConnectivityList)                
        Coords = cartesianToBarycentric(TR,TRGidx,PT);
        if min(Coords) >= 0
            W(NewInterps(PTidx),TR(TRGidx,:)) = Coords;%THIRD ARGUEMENT 
            break
        end   
    end
end

%Interpolate
Fun  =  W*K;