function [v] = ParalleTransportVectorField(S,TransMtx,EdgeInf,seedtrg,seedV)

num_trg = size(S.trg,1);
init_trg = seedtrg;
v = zeros(3,num_trg); %%%%%
v(:,init_trg) = seedV;

indicator_trg = zeros(num_trg,1);
indicator_trg(init_trg) = 1;
Red = [init_trg];  %%%%%  New added triangle
Green = [1:num_trg];  %%%%%  Triangles have been touched
Black = [init_trg];  %%%%% triangles have been touched
newRed = [];

while length(Black) < num_trg
    
    for i = 1:length(Red)
        f = Red(i);
        p = [S.trg(f,:) S.trg(f,1)];  %%%% vertex indices of triangle f
        for k = 1:3
            count = EdgeInf.edgeidx(p(k),p(k+1));
            f1 = EdgeInf.edgeadjface(count,1);
            f2 = EdgeInf.edgeadjface(count,2);
            if f1 == 0
                f1 = f2;
            elseif f2 == 0;
                f2 = f1;
            end
            if (f1 == f && indicator_trg(f2) == 0)
                v(:,f2) = TransMtx{count}*v(:,f);
                indicator_trg(f2) = 1;
                newRed = union(newRed,f2);
            elseif (f2 == f && indicator_trg(f1) == 0)
                v(:,f1) = TransMtx{count}'*v(:,f);
                indicator_trg(f1) = 1;
                newRed = union(newRed,f1);
            end
        end
        %Red = setdiff(Red,f);
        Black = union(Black,f);
    end
    Red = newRed;
    newRed = [];
    
end
