function [ wDenom ] = updateCache( sp )
    Temp=sp.OtfProvider.otf;
    siz=size(Temp(:,:,1));
    w=siz(2);
    h=siz(1);
    wDenom=zeros(2*h,2*w);
    for d=1:3
        for b=1:2
            wd=wDenom;
            [wDenom]=addWienerDenominator_2D(wd,sp,d,b);
        end
    end
end