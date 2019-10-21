function [Surf] = DefineMassAndStiff(Surf,eigs)

if nargin == 1
    eigs = 0;
end

[V, D, MV, MM, SM] = LBOeigs(Surf,eigs);

if eigs == 1
    Surf.V = V;
    Surf.D = D;
end
Surf.M      = MM;
Surf.S      = SM;
Surf.Mdiga  = MV;
Surf.nPts   = length(Surf.pts);
Surf.nTrg   = length(Surf.trg);