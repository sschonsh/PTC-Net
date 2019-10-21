clear

%% Data
fname = 'mesh018';
opts.print = 1;
[surf.pt, surf.trg] = ReadOFF([fname '.off']);
surf.pts = surf.pt;
num_pt = size(surf.pt,1);
num_trg = size(surf.trg,1);
scale = 2;
Seedpt = 8095; 
pt2    = 9567;


% z = 3*cos(6*surf.pt(:,1)).*sin(8*surf.pt(:,2));
% surf.pt(:,3) = z;
% surf.pt(:,1:2) = surf.pt(:,1:2)*10;

%KNN stuff
Locals = knnsearch(surf.pt,surf.pt,'K',350);

%% Define Vector Field and frame

%on Triangle
TrgNormal = findUpNormal(surf);
for i = 1:length(surf.trg)
    %find bary and area
    p1 = surf.pt(surf.trg(i,1),:);
    p2 = surf.pt(surf.trg(i,2),:);
    p3 = surf.pt(surf.trg(i,3),:);
    v1 = p1-p2;
    v2 = p1-p3;
    surf.area(i)   = .5*norm(cross(v1,v2));
    surf.bary(i,:) = (p1+p2+p3)/3;
end

%set Primary V
PrimaryV = [1,0,0];

%caluate normals as weighted average from local faces
hold on;
Normals = zeros(length(surf.pt),3);
V1 = Normals;
V2 = Normals;
for i = 200:length(surf.pt)%2--
    [FirstRingTRG,~] = find(surf.trg == i);
    weights      = bsxfun(@times,surf.area(FirstRingTRG),ones(3,length(FirstRingTRG)))'; 
    WNorms       = TrgNormal(FirstRingTRG,:).*weights;
    TempNormal   = sum(WNorms);
    Normals(i,:) = TempNormal/norm(TempNormal);
    
    %Define Vs
    TempV1  = PrimaryV -(PrimaryV*Normals(i,:)')*Normals(i,:);
    V1(i,:) = TempV1/norm(TempV1);
    V2(i,:) = cross(Normals(i,:),V1(i,:));

end

%% Project onto tangent plane at point
p   = surf.pt(Seedpt,:);
TN  = Normals(Seedpt,:);
figure
hold on
for i = 1:length(Locals(Seedpt,:))
    %project onto plane
    q = surf.pt(Locals(Seedpt,i),:);
    proj = q - ((q-p)*TN')*TN;
    v = proj - p;
    
    %record coeffiects
    Coeff1(i) = v*V1(Seedpt,:)';
    Coeff2(i) = v*V2(Seedpt,:)';
    Coeff3(i) = v*TN';%should be zeros
    
end

%% Define Kernal
K1 = zeros(length(surf.pt),1);
K1(Locals(Seedpt,:)) = (15./exp(400*(Coeff1.^2+Coeff2.^2))).*cos(atan(Coeff2./(Coeff1+eps)))-1.9330e-05;
ViewMesh(surf,K1)
colorbar

%% Iteropolate

%2
InterpFun2 = scatteredInterpolant(Coeff1',Coeff2',K1(Locals(Seedpt,:)));
p   = surf.pt(pt2,:);
TN  = Normals(pt2,:);
for k = 1:length(Locals(pt2,:))
    %project onto plane
    q = surf.pt(Locals(pt2,k),:);
    proj = q - ((q-p)*TN')*TN;
    v = proj - p;
    LCoeff1(k) = v*V1(pt2,:)';
    LCoeff2(k) = v*V2(pt2,:)';
end
Values = InterpFun2(LCoeff1',LCoeff2');
K2 = zeros(length(surf.pt),1);
K2(Locals(pt2,:)) = Values;

%3
InterpFun3 = scatteredInterpolant(2*Coeff1',2*Coeff2',K1(Locals(Seedpt,:)));
pt2 = Seedpt;
p   = surf.pt(pt2,:);
TN  = Normals(pt2,:);
for k = 1:length(Locals(pt2,:))
    %project onto plane
    q = surf.pt(Locals(pt2,k),:);
    proj = q - ((q-p)*TN')*TN;
    v = proj - p;
    LCoeff1(k) = v*V1(pt2,:)';
    LCoeff2(k) = v*V2(pt2,:)';
end
Values = InterpFun3(LCoeff1',LCoeff2');
K3 = zeros(length(surf.pt),1);
K3(Locals(pt2,:)) = Values;


%4
Ang = pi/2;
InterpFun4 = scatteredInterpolant(Coeff1',Coeff2',K1(Locals(Seedpt,:)));
pt2 = Seedpt;
p   = surf.pt(pt2,:);
TN  = Normals(pt2,:);
for k = 1:length(Locals(pt2,:))
    %project onto plane
    q = surf.pt(Locals(pt2,k),:);
    proj = q - ((q-p)*TN')*TN;
    v = proj - p;
    
    %find rotated frame frame
    CofB = [V1(pt2,:); V2(pt2,:); TN];
    Rot  = [cos(Ang), -sin(Ang), 0;...
           sin(Ang),  cos(Ang), 0;...
           0       , 0        , 0;];
    B1 = CofB'*Rot*[1,0,0]'; %CofB'*Rot*CofB*(V1(i,:))';
    B2 = CofB'*Rot*[0,1,0]'; %CofB'*Rot*CofB*(V2(i,:))';
    
    LCoeff1(k) = v*B1;
    LCoeff2(k) = v*B2;
end
Values = InterpFun4(LCoeff1',LCoeff2');
K4 = zeros(length(surf.pt),1);
K4(Locals(pt2,:)) = Values;



%5
pt2 = 7217;
Ang = pi/2;
InterpFun5 = scatteredInterpolant(2*Coeff1',2*Coeff2',K1(Locals(Seedpt,:)));
p   = surf.pt(pt2,:);
TN  = Normals(pt2,:);
for k = 1:length(Locals(pt2,:))
    %project onto plane
    q = surf.pt(Locals(pt2,k),:);
    proj = q - ((q-p)*TN')*TN;
    v = proj - p;
    
    %find rotated frame frame
    CofB = [V1(pt2,:); V2(pt2,:); TN];
    Rot  = [cos(Ang), -sin(Ang), 0;...
           sin(Ang),  cos(Ang), 0;...
           0       , 0        , 0;];
    B1 = CofB'*Rot*[1,0,0]'; %CofB'*Rot*CofB*(V1(i,:))';
    B2 = CofB'*Rot*[0,1,0]'; %CofB'*Rot*CofB*(V2(i,:))';
    
    LCoeff1(k) = v*B1;
    LCoeff2(k) = v*B2;
end
Values = InterpFun5(LCoeff1',LCoeff2');
K5 = zeros(length(surf.pt),1);
K5(Locals(pt2,:)) = Values;



%% Plots
figure
subplot(1,5,1)
ViewMesh(surf,K1)
xlabel('X','FontSize',12)
ylabel('Y','FontSize',12)
%axis on
title('Original Kernal','FontSize',18)

subplot(1,5,2)
ViewMesh(surf,K2)
xlabel('X','FontSize',12)
ylabel('Y','FontSize',12)
title('Transltion','FontSize',18)
%axis on

subplot(1,5,3)
ViewMesh(surf,K3)
xlabel('X','FontSize',12)
ylabel('Y','FontSize',12)
title('Dilation','FontSize',18)
%axis on

subplot(1,5,4)
ViewMesh(surf,K4)
xlabel('X','FontSize',12)
ylabel('Y','FontSize',12)
title('Rotation','FontSize',18)
%axis on

subplot(1,5,5)
ViewMesh(surf,K5)
xlabel('X','FontSize',12)
ylabel('Y','FontSize',12)
title('All three','FontSize',18)
%axis on

%% MOVIE

B = surf.trg;
G = digraph([B(:,1);B(:,2);B(:,3)],[B(:,2);B(:,3);B(:,1)]);
A = adjacency(G);


%OPTIONS

% PointIdx = [shortestpath(G,8095,9567),shortestpath(G,9567,8095)];
% AngIdx   = zeros(size(PointIdx));
% ScaleIdx    = ones(size(PointIdx));

% PointIdx = 8095*ones([100,1]);
% AngIdx   = linspace(0,3*pi,100);
% ScaleIdx    = ones(size(PointIdx));

PointIdx = [shortestpath(G,8095,9567),shortestpath(G,9567,8095)];
AngIdx   = linspace(0,4*pi,length(PointIdx));
ScaleIdx = [linspace(1,.4,floor(length(PointIdx)/4)), linspace(.4,1.3,floor(length(PointIdx)/2)),linspace(1.3,1,ceil(length(PointIdx)/4))];


%Big Loop over frames
Frames  = zeros(length(PointIdx),length(surf.pt));
figure
for i = 1:length(PointIdx)
    pt2 = PointIdx(i);
    Ang = AngIdx(i);
    Scale = ScaleIdx(i);
    p   = surf.pt(pt2,:);
    TN  = Normals(pt2,:);

    for k = 1:length(Locals(pt2,:))
        %project onto plane
        q = surf.pt(Locals(pt2,k),:);
        proj = q - ((q-p)*TN')*TN;
        v = proj - p;

        %find rotated frame frame
        CofB = [V1(pt2,:); V2(pt2,:); TN];
        Rot  = [cos(Ang), -sin(Ang), 0;...
               sin(Ang),  cos(Ang), 0;...
               0       , 0        , 0;];
        B1 = CofB'*Rot*[1,0,0]'; %CofB'*Rot*CofB*(V1(i,:))';
        B2 = CofB'*Rot*[0,1,0]'; %CofB'*Rot*CofB*(V2(i,:))';

        LCoeff1(k) = v*B1;
        LCoeff2(k) = v*B2;
    end
    Values = InterpFun4(Scale*LCoeff1',Scale*LCoeff2');
    K = zeros(length(surf.pt),1);
    K(Locals(pt2,:)) = Values;
    
    %Plot
    clf
    ViewMesh(surf,K)
    Frames(i,:) = K;
    pause(.1)


    
end









