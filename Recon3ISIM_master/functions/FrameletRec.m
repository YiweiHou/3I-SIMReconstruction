function Rec = FrameletRec(C,L)
    R{1} = [1 2 1]/4;
    R{2} = [-1 0 1]/4*sqrt(2);
    R{3} = [-1 2 -1]/4;
for k = L:-1:2
    C{k-1}{1,1} = FraRec(C{k},R,k);
end
Rec = FraRec(C{1},R,1);
end

function Rec = FraRec(C,R,L)
nR = length(R);
ImSize = size(C{1,1});
Rec = zeros(ImSize);

for i = 1:nR
    temp = zeros(ImSize);
    for j = 1:nR
        M2 = R{j};
        temp = temp + (ConvSymAsym((C{i,j})',M2,L))';
    end
    M1 = R{i};
    Rec = Rec + ConvSymAsym(temp,M1,L);
end
end