function Surf = DefineDownSample(Surf,nPts,seed)

DownPts(1) = seed;
close
figure;
ViewMesh(Surf.pts,Surf.trg);
hold on
plot3(Surf.pts(DownPts(1),1),Surf.pts(DownPts(1),2),Surf.pts(DownPts(1),3),'r*','LineWidth',3)
hold on

% compute new longest distance
Distances = Surf.Distances;
for j = 1:nPts-1    
     minDistances = min(Distances(DownPts,:),[],1);
    [~,IDX] = max(minDistances);
    DownPts(j+1) = IDX;

%      plot3(Surf.pts(DownPts(j+1),1),Surf.pts(DownPts(j+1),2),Surf.pts(DownPts(j+1),3),'r*','LineWidth',3)
%      pause(0)
end

Surf.DownSample = sort(DownPts);