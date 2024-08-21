function [para] = Rounding(para)
if para>=1e-2
    count=0;
    while para<1
        para=para*10;
        count=count+1;
    end
    para=para*10;
    para=round(para);
    para=para/10^(count+1);
else
        count=0;
    while para<1
        para=para*10;
        count=count+1;
    end
    para=round(para);
    para=para/10^(count);
end
end

