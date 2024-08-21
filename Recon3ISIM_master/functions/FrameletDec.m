function Dec = FrameletDec(A,L)
    D{1} = [1 2 1]/4;
    D{2} = [1 0 -1]/4*sqrt(2);
    D{3} = [-1 2 -1]/4;
    kDec = A;
for k = 1:L
    Dec{k} = FraDec(kDec,D,k);
    kDec = Dec{k}{1,1};
end
end

function Dec = FraDec(A,D,L)
nD = length(D);
for i = 1:nD
    M1 = D{i};
    tempi = ConvSymAsym(A,M1,L);
    for j = 1:nD
        M2 = D{j};        
        tempj = ConvSymAsym(tempi',M2,L);
        Dec{i,j} = single(tempj');
    end
end
end
