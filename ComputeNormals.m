function [Surf] = ComputeNormals(Surf);

pts = Surf.pts;
trg = Surf.trg;

TrgNormal = findUpNormal(Surf);
for i = 1:length(Surf.trg)
    %find bary and area
    p1 = pts(trg(i,1),:);
    p2 = pts(trg(i,2),:);
    p3 = pts(trg(i,3),:);
    v1 = p1-p2;
    v2 = p1-p3;
    Surf.area(i,1)   = .5*norm(cross(v1,v2));
    Surf.bary(i,:) = (p1+p2+p3)/3;
end

%caluate normals as weighted average from local faces
Normals = zeros(length(Surf.pts),3);
for i = 1:length(Surf.pts)
    [FirstRingTRG,~] = find(Surf.trg == i);
    weights      = bsxfun(@times,Surf.area(FirstRingTRG)',ones(3,length(FirstRingTRG)))'; 
    WNorms       = TrgNormal(FirstRingTRG,:).*weights;
    TempNormal   = sum(WNorms);
    Normals(i,:) = TempNormal/norm(TempNormal);
end
Surf.normals = Normals;
