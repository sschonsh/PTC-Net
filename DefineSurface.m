SeedPts = [1,100,222,132];


%Create Surface
N = 28;
Surf.nPts = N^2;
[x,y]    = meshgrid(linspace(0,10,N),linspace(0,10,N));
pts(:,1) = reshape(x,N^2,1);
pts(:,2) = reshape(y,N^2,1);
pts(:,3) = 2*cos(.6*pts(:,1)).*sin(.8*pts(:,2));
Surf.pts = pts;
Surf.trg = delaunay(pts(:,1:2));

%VFs
Surf = DefineMassAndStiff(Surf,1);
Surf = ComputeNormals(Surf);
Surf = ComputeGeodesicVFs(Surf,SeedPts);
Surf = ComputeCanonicalVFs(Surf,[1,1,1]);

%Stencil
Radius  = 3;
NumRing = 2;
Surf = SimpleStencil(Surf,Radius,NumRing);


%Interp
Surf = DefineLocalInfoWithRestriction(Surf);
Surf = StackRandWNoPad(Surf);

%Save
sname = ['ExSurf.mat'];
save(sname,'Surf') 
