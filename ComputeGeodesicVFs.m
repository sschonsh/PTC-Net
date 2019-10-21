function [Surf] = ComputeGeodesicVFs(Surf,Seedpts)

%Define number of VFs to make
Surf.nVFs = length(Seedpts);

for i = 1:length(Seedpts)
    
    %% define vector fields
    opts.print = 0;
    VecSeedpt1 = Seedpts(i);
    opts.init  = VecSeedpt1;
    [TransMtx1, LocalFrame1, EdgeInf1, Surf] = ConstructTransMtx(Surf,opts);
    a = [1,1];
    Seedtrg = 1;
    SeedV = a(1)*LocalFrame1.v1(:,Seedtrg) + a(2)*LocalFrame1.v1(:,Seedtrg);  %%%%%%% 3D column vector
    TrgV = ParalleTransportVectorField(Surf,TransMtx1,EdgeInf1,Seedtrg,SeedV)';
    
    for j = 1:Surf.nPts
        [FirstRingTRG,~] = find(Surf.trg == j);
        weights = bsxfun(@times,Surf.area(FirstRingTRG)',ones(3,length(FirstRingTRG)))'; 
        WVs     = TrgV(FirstRingTRG,:).*weights;
        TempV   = sum(WVs,1);
        V(j,:)  = TempV/norm(TempV);
    end
    Surf.VF{i} = V;
end