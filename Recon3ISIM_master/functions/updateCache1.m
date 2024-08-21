function [ wDenom ] = updateCache1( sp )
    wDenom=zeros(2*sp.imgSize,2*sp.imgSize);
    for d=1:3
        wDenom=wDenom+1.0*WienerDenominatorNoAtt_2D(sp,d,1)+1.0*WienerDenominatorUseAtt_2D(sp,d,2);
    end
end