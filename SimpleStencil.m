function [Surf] = SimpleStencil(Surf,MaxRadius,NumRing)

Points  = [0,0];
PerRing  = 4;
Radius   = 0;
RingInc  = MaxRadius/NumRing;
for i = 1:NumRing
   Radius = Radius + RingInc;
   Theta  = linspace(0,2*pi,PerRing+1);
   for j = 1:PerRing
       NewPoint = Radius*[sin(Theta(j)),cos(Theta(j))];
       Points = [Points; NewPoint];
   end
   PerRing = 2*PerRing;
end
Surf.nTempPts = length(Points);
Surf.Template = Points;
Surf.TempRad  = MaxRadius;
