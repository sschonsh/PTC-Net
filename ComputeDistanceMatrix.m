function [Surf] = ComputeDistanceMatrix(Surf)

Distances = zeros(length(Surf.pts));

for i = 1:length(Surf.pts)
    [D,~,~]  = perform_fast_marching_mesh(Surf.pts, Surf.trg, i);
    Distances(:,i) = D;
end

Surf.Distances = Distances;