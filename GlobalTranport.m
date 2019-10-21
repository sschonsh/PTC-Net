function M = GlobalTranport(K,Surf,VFidx)

M = sparse(Surf.nPts,Surf.nPts);
for i = 1:Surf.nPts
   Temp   =  Surf.W{VFidx,i}*K;
   M(:,i) =  Temp;
end
    