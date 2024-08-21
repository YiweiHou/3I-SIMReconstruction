function C = ConvSymAsym(A,M,L)
[m,n] = size(A);
nM = length(M);
step = 2^(L-1);
ker = zeros(step*(nM-1)+1,1);
ker(1:step:step*(nM-1)+1,1) = M;
lker = floor(length(ker)/2);
C = imfilter(A,ker,'circular');
end