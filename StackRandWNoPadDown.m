function [Surf] = StackRandWNoPadDown(Surf,IDX)


for vf = 1:Surf.nVFs
Wfull = [];
Rfull = [];
Ifull = [];
count = 1;
%start = 1;
    for i = 1:length(IDX)
        pt = IDX(i);
        Wcur = Surf.W{vf,pt};
        Wfull = [Wfull; Wcur];

        Rcur = Surf.R{pt};
        Rfull = [Rfull, Rcur];

        [I,~] = find(Rcur);
        Surf.I{pt} = I;
        Surf.Isize{pt} = length(I);
%         for j = 1:length(I)
%             IRowIndex(count) = pt;
%             count = count+1;
%         end
         Ifull = [Ifull; I];
%         [A,~] = size(Wfull);
    end
    %Save
    if vf == 1
        Surf.Wfull1 = sparse(Wfull);
        Surf.Rfull1 = sparse(Rfull);
        Surf.Ifull1 = Ifull;
        [Surf.IRowIndex1, Surf.IColumnIndex1] = find(Surf.Rfull1);
%         Surf.IRowIndex1    = IRowIndex;
%         Surf.IColumnIndex1 = 1:length(Ifull);
        Surf.L1 = length(Ifull);
    elseif vf == 2
        Surf.Wfull2 = sparse(Wfull);
        Surf.Rfull2 = sparse(Rfull);
        Surf.Ifull2 = Ifull;
        [Surf.IRowIndex2, Surf.IColumnIndex2] = find(Surf.Rfull2);
        Surf.L2 = length(Ifull);
    elseif vf == 3
        Surf.Wfull3 = sparse(Wfull);
        Surf.Rfull3 = sparse(Rfull);
        Surf.Ifull3 = Ifull;
        [Surf.IRowIndex3, Surf.IColumnIndex3] = find(Surf.Rfull3);
        Surf.L3 = length(Ifull);
    elseif vf == 4
        Surf.Wfull4 = sparse(Wfull);
        Surf.Rfull4 = sparse(Rfull);
        Surf.Ifull4 = Ifull;
        [Surf.IRowIndex4, Surf.IColumnIndex4] = find(Surf.Rfull4);
        Surf.L4 = length(Ifull);
    else
        disp('Error, only support up to 4vf right now')
    end
end