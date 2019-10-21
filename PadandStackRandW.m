function [Surf] = PadandStackRandW(Surf)

L = 0;
for i = 1:Surf.nPts
   [a,~] = size(Surf.W{i});
   if a > L
       L = a;
   end
end

Surf.maxL = L;


for j = 1:Surf.nVFs
    Wfull = [];
    Rfull = [];
    for i = 1:Surf.nPts
        [a,~] = size(Surf.W{j,i});
        Wpad = Surf.W{j,i};
        Wpad(a+1:L,:) = 0; %= padarray(Surf.W{i},L-b,0);
        Wfull = [Wfull; Wpad];
    
        [~,b] = size(Surf.R{i});
        Rpad = Surf.R{i};
        Rpad(:,b+1:L) = 0;
        Rfull = [Rfull, Rpad]; 
    end


Surf.Wfull{j} = sparse(Wfull);
Surf.Rfull{j} = sparse(Rfull);

end