%This demos how to set up a surf with everything you need to define a
%network
clear
%%%%%THIS NEEDS TO BE RUN IN COMMAND WINDOW ONCE TO COMPILE MEX FILES
% cd toolbox_fast_marching
% compile_mex 
% cd ..

%% Load Data
Surf = DefineSimpleManifold;
subplot(1,3,1);
ViewMesh(Surf.pts,Surf.trg);
title('Mesh')

%% Define Intrinsics
Surf = DefineMassAndStiff(Surf);

%% Define VFs
Surf = ComputeNormals(Surf);
%Use Geodesic
SeedPts = [1,28^2];%Seed points for geodesic vector Fields
Surf = ComputeGeodesicVFs(Surf,SeedPts);


%% Choose Template for Conv
Radius  = 1/9;
NumRing = 3;
Surf = SimpleStencil(Surf,Radius,NumRing);
subplot(1,3,2)
plot(Surf.Template(:,1),Surf.Template(:,2),'r+','LineWidth',3);
title('Stencil')

%% Define Neighborhoods and Coefficents and Weights
Surf = DefineLocalInfo(Surf);

%% Create a Random Kernel
K = randn(length(Surf.Template),1);

%% Compute Global Transport Matrix
VFidx = 1;%Choose which of the vecotr fields to use
M = GlobalTranport(K,Surf,VFidx);

%% Show some tranportaions
subplot(1,3,3)
ViewMesh(Surf.pts,Surf.trg,M(:,65)+M(:,366)+M(:,end)+M(:,501));
hold on;
plot3(Surf.pts(65 ,1),Surf.pts(65 ,2),Surf.pts(65 ,3),'k*','MarkerSize',8)
plot3(Surf.pts(366,1),Surf.pts(366,2),Surf.pts(366,3),'k*','MarkerSize',8)
plot3(Surf.pts(end,1),Surf.pts(end,2),Surf.pts(end,3),'k*','MarkerSize',8)
plot3(Surf.pts(501,1),Surf.pts(501,2),Surf.pts(501,3),'k*','MarkerSize',8)
title('Transportations')

%% Tranport to spot with scale + rotation
location1 = 200;
location2 = 555;
scale    = 2;
rotation = pi;
Kprime  = SpotTranport(K,Surf,VFidx,location1,scale,rotation);
Kprime2 = SpotTranport(K,Surf,VFidx,location1,1,0);
Kprime3 = SpotTranport(K,Surf,VFidx,location2,1,0);
Kprime4 = SpotTranport(K,Surf,VFidx,location2,.75,0);

figure;
subplot(2,2,1)
ViewMesh(Surf.pts,Surf.trg,Kprime)
hold on
plot3(Surf.pts(location1,1),Surf.pts(location1,2),Surf.pts(location1,3),'k*','MarkerSize',10)
title(['Location:' num2str(location1) ', Scale:2'])

subplot(2,2,2)
ViewMesh(Surf.pts,Surf.trg,Kprime2)
hold on
plot3(Surf.pts(location1,1),Surf.pts(location1,2),Surf.pts(location1,3),'k*','MarkerSize',10)
title(['Location:' num2str(location1) ', Scale:1'])

subplot(2,2,3)
ViewMesh(Surf.pts,Surf.trg,Kprime3)
hold on
plot3(Surf.pts(location2,1),Surf.pts(location2,2),Surf.pts(location2,3),'k*','MarkerSize',10)
title(['Location:' num2str(location2) ', Scale:1'])

subplot(2,2,4)
ViewMesh(Surf.pts,Surf.trg,Kprime4)
hold on
plot3(Surf.pts(location2,1),Surf.pts(location2,2),Surf.pts(location1,3),'k*','MarkerSize',10)
title(['Location:' num2str(location2) ', Scale:0.75'])