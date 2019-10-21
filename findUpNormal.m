function N = findUpNormal(surf)

pts = surf.pts;
trg = surf.trg;
N = zeros(size(trg));

for i = 1:length(trg)
    p1 = pts(trg(i,1),:);
    p2 = pts(trg(i,2),:);
    p3 = pts(trg(i,3),:);
    
    v1 = p1-p2;
    v2 = p1-p3;
     
    ntemp = cross(v1,v2);
    ntemp = ntemp./(norm(ntemp));
    
    if ntemp(3) < 0
        ntemp = -ntemp;
    end
    
    N(i,:) = ntemp;
    
end