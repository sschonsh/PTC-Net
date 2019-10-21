function [TransMtx, LocalFrame, EdgeInf, Surf] = ConstructTransMtx(Surf,opts)

%%
pt = Surf.pts';
trg = Surf.trg;

num_pt = length(pt);
num_trg = size(trg,1);
edges = unique([trg(:,1:2);trg(:,2:3)],'rows');%
num_edge = length(edges);%
%num_edge = 3*num_trg/2;
facenormal = zeros(3,num_trg);


EdgeInf.rowidx = zeros(num_edge,1);
EdgeInf.colidx = zeros(num_edge,1);
EdgeInf.edgeidx = sparse(num_pt,num_pt);
EdgeInf.edgeadjface = zeros(num_edge,2);
edge = sparse(num_pt,num_pt);
ptaroundpt = cell(num_pt,1);
ptaroundtrg = cell(num_pt,1);
LocalFrame.v1 = zeros(3,num_trg);  %%%% per triangle
LocalFrame.v2 = zeros(3,num_trg);
trgcenter = zeros(3,num_trg);
count = 0;
for i = 1:num_trg
    p1 = trg(i,1);
    p2 = trg(i,2);
    p3 = trg(i,3);
    %%EdgeVect(i,:) = [norm(pt(p2,:) - pt(p1,:))^2, norm(pt(p3,:) - pt(p1,:))^2,norm(pt(p2,:) - pt(p3,:))^2];
    v1 = pt(:,p1);   v2 = pt(:,p2);   v3 = pt(:,p3);
    trgcenter(:,i) = (v1 + v2 + v3)/3;
    tempnormal = cross(v2-v1,v3-v1);
    tempnormal = tempnormal/norm(tempnormal);
    facenormal(:,i) = tempnormal;
    ptaroundpt{p1} = union(ptaroundpt{p1},[p2 p3]);
    ptaroundpt{p2} = union(ptaroundpt{p2},[p1 p3]);
    ptaroundpt{p3} = union(ptaroundpt{p3},[p1 p2]);
    ptaroundtrg{p1}=union(ptaroundtrg{p1},i);
    ptaroundtrg{p2}=union(ptaroundtrg{p2},i);
    ptaroundtrg{p3}=union(ptaroundtrg{p3},i);
    
    LocalFrame.v1(:,i) = (v2 - v1)/norm(v2 - v1);
    LocalFrame.v2(:,i) = cross(tempnormal,LocalFrame.v1(:,i));
    
    
    if edge(p1,p2) == 0
        count = count + 1;
        EdgeInf.rowidx(count) = p1;
        EdgeInf.colidx(count) = p2;
        EdgeInf.edgeidx(p1,p2) = count;
        EdgeInf.edgeidx(p2,p1) = count;
        EdgeInf.edgeadjface(count,1) = i;
        edge(p1,p2) = i; edge(p2,p1) = -1;
    else
        EdgeInf.edgeadjface(EdgeInf.edgeidx(p1,p2),2) = i;
        edge(p1,p2) = i;
    end
    
    if edge(p2,p3) == 0
        count = count + 1;
        EdgeInf.rowidx(count) = p2;
        EdgeInf.colidx(count) = p3;
        EdgeInf.edgeidx(p2,p3) = count;
        EdgeInf.edgeidx(p3,p2) = count;
        EdgeInf.edgeadjface(count,1) = i;
        edge(p2,p3) = i; edge(p3,p2) = -1;
    else
        EdgeInf.edgeadjface(EdgeInf.edgeidx(p2,p3),2) = i;
        edge(p2,p3) = i;
    end
    
    if edge(p3,p1) == 0
        count = count + 1;
        EdgeInf.rowidx(count) = p3;
        EdgeInf.colidx(count) = p1;
        EdgeInf.edgeidx(p3,p1) = count;
        EdgeInf.edgeidx(p1,p3) = count;
        EdgeInf.edgeadjface(count,1) = i;
        edge(p3,p1) = i; edge(p1,p3) = -1;
    else
        EdgeInf.edgeadjface(EdgeInf.edgeidx(p3,p1),2) = i;
        edge(p3,p1) = i;
    end
    
end

EdgeInf.num_edge = count;
Surf.facenormal = facenormal;
Surf.aroundtrg = ptaroundtrg;
Surf.trgcenter = trgcenter;

ptaroundtrgpath = cell(num_pt,1);
aroundedge = cell(num_pt,1);
aroundpt = cell(num_pt,1);
for i = 1:num_pt
    temparoundtrg = [full(edge(ptaroundpt{i},i)),full(edge(i,ptaroundpt{i}))'];
    temppath = temparoundtrg(1,:);
    tempedge = [i ptaroundpt{i}(1)];
    while temppath(end) ~= temppath(1)
        tempidx = find(temparoundtrg(:,1) == temppath(end));
        addtrgidx = temparoundtrg(tempidx,2);
        temppath = [temppath addtrgidx];
        if  ~isempty(ptaroundpt{i}(tempidx))%
            tempedge = [tempedge; [i ptaroundpt{i}(tempidx)]];%
        end%
    end
    ptaroundtrgpath{i} = temppath;
    aroundedge{i} = tempedge;
    aroundpt{i} = tempedge(:,2);
end
%
%
Surf.ptaroundtrgpath = ptaroundtrgpath;
Surf.aroundedge = aroundedge;
Surf.aroundpt = aroundpt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Fast Marching for distance
initialpt = opts.init;
D = FastMarching(initialpt,Surf);
%%%%%  Construct Local frame

for i = 1:num_trg
    E21 = pt(:,trg(i,2)) - pt(:,trg(i,1));
    E31 = pt(:,trg(i,3)) - pt(:,trg(i,1));
    d21 = D(trg(i,2)) - D(trg(i,1));
    d31 = D(trg(i,3)) - D(trg(i,1));
    
    shapegrad = [E21,E31]*([E21'*E21,E21'*E31;E21'*E31,E31'*E31]\[d21;d31]);
    LocalFrame.v1(:,i) = shapegrad/norm(shapegrad);
    LocalFrame.v2(:,i) = cross(facenormal(:,i),LocalFrame.v1(:,i));
end

%%%%%%%%   Transition matrix between two adjacency triangles
%%
num_edge = EdgeInf.num_edge;
TransMtx = cell(num_edge,1);

for i = 1:num_edge
    p1 = EdgeInf.rowidx(i);
    p2 = EdgeInf.colidx(i);
    e = pt(:,p2) - pt(:,p1); e = e/norm(e);
    count = EdgeInf.edgeidx(p1,p2);
    f1 = EdgeInf.edgeadjface(count,1);
    f2 = EdgeInf.edgeadjface(count,2);
    
    
    %%%%%%%  use rotation
    if f1 ~= 0 && f2 ~= 0;
        T_f1 = Rotation(LocalFrame.v1(:,f1),e,facenormal(:,f1)); %%% Transform f1.v1 to e
        T_f2 = Rotation(e, LocalFrame.v1(:,f2),facenormal(:,f2)); %%% Transform e to f2.v1
        T_edge = Rotation(facenormal(:,f1),facenormal(:,f2),e);  %%%% Transform f1.n to f2.n
    elseif f1 == 0;
        T_f1 = Rotation(LocalFrame.v1(:,f2),e,facenormal(:,f2)); %%% Transform f1.v1 to e
        T_f2 = Rotation(e, LocalFrame.v1(:,f2),facenormal(:,f2)); %%% Transform e to f2.v1
        T_edge = Rotation(facenormal(:,f2),facenormal(:,f2),e);  %%%% Transform f1.n to f2.n     
    elseif f2 == 0;
        T_f1 = Rotation(LocalFrame.v1(:,f1),e,facenormal(:,f1)); %%% Transform f1.v1 to e
        T_f2 = Rotation(e, LocalFrame.v1(:,f1),facenormal(:,f1)); %%% Transform e to f2.v1
        T_edge = Rotation(facenormal(:,f1),facenormal(:,f1),e);  %%%% Transform f1.n to f2.n
    end
    
    TransMtx{i} = T_f2*T_edge*T_f1;  %%% TransMtx from f1 to f2
    
end
