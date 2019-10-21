function [Surf] = ComputeCanonicalVFs(Surf,Vec)

%Define number of VFs to make
[Surf.nVFs,~] = size(Vec);

for i = 1:Surf.nVFs
    PrimaryV = Vec(i,:);
    %% define vector fields
     for j = 1:length(Surf.pts)         
        TempV1  = PrimaryV -(PrimaryV*Surf.normals(j,:)')*Surf.normals(j,:);
        V(j,:) = TempV1/norm(TempV1);
    end
    Surf.VF{i} = V;
end