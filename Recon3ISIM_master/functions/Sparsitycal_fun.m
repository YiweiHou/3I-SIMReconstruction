function [sparsity_framelet,sparsity_curvelet] = Sparsitycal_fun(f)
f=single(f);
f=f/max(max(f));
    wav1=@(U)FrameletDec(U,2);
    iwav1=@(C)FrameletRec(C,2);
    wav2 = @(X)CurveletDec(X,1,2); 
    iwav2 =@(C)real(CurveletRec(C,1));

FrameletC=wav1(f);
CurveletC=wav2(f);

FrameletC=FrameletC{1,1};
[m,n]=size(FrameletC);
Frameletco=[];

%% Framelet Gini index
total=0;
for i=1:m
    for j=1:n
        tmp=FrameletC{i,j};
        tmp=abs(tmp);
        total=total+sum(sum(tmp));
        Frameletco=[Frameletco,tmp(:)'];
    end
end
Frameletco=sort(Frameletco);
[~,N]=size(Frameletco);
seq=1:1:N;
sequence=(N-seq+1/2)/N/total;
term1=2*Frameletco.*sequence;
term1=sum(term1);
sparsity_framelet=1-term1;


%% Curvelet Gini index
total=0;
    [m,n]=size(CurveletC);
    Curveletco=[];
    for ii=1:m
        for jj=1:n
            tmp1=CurveletC{ii,jj};
            z=length(tmp1);
            for kk=1:z
                tmp=tmp1{1,kk};
                tmp=abs(tmp);
                total=total+sum(sum(tmp));
                Curveletco=[Curveletco,tmp(:)'];
            end
        end
    end
    Curveletco=sort(Curveletco);
[~,N]=size(Curveletco);
seq=1:1:N;
sequence=(N-seq+1/2)/N/total;
term1=2*Curveletco.*sequence;
term1=sum(term1);
sparsity_curvelet=1-term1;
end

