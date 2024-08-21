function [ret] = fromEstimate(ret)
    ret.isMultiband=0;
    ret.isEstimate=1;
    vals1=zeros(1,ret.sampleLateral);
    valsAtt=zeros(1,ret.sampleLateral);
    valsOnlyAtt=zeros(1,ret.sampleLateral);


    for I=1:ret.sampleLateral
        v=abs(I-1)/ret.sampleLateral;
        r1=valIdealOTF(v)*power(ret.estimateAValue,v);
        vals1(I)=r1;
    end

    for I=1:ret.sampleLateral
        dist=abs(I-1)*ret.cyclesPerMicron;
        valsOnlyAtt(I)=valAttenuation(dist,ret.attStrength,ret.attFWHM);
        valsAtt(I)=vals1(I)*valsOnlyAtt(I);
    end

    ret.vals=vals1;

    ret.valsAtt=valsAtt;
    ret.valsOnlyAtt=valsOnlyAtt;

end