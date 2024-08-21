function [ ret ] = Airy_giveOtf(param,a)
    ret.na = param.NA;
    ret.lambda = param.lambda;
    ret.cutoff = 1000/(0.5*ret.lambda/ret.na);
    
    ret.imgSize = param.imgSize;
    ret.cyclesPerMicron = param.cyclesPerMicron;
    ret.sampleLateral = ceil(ret.cutoff/ret.cyclesPerMicron)+1;

    ret.estimateAValue = a;
    ret.maxBand = 2;
    ret.attStrength = param.attStrength;
    ret.attFWHM = 1.0;    
    ret.useAttenuation = 1;

    ret=fromEstimate(ret);
    
    ret.otf=zeros(param.imgSize,param.imgSize);
    ret.otfatt=zeros(param.imgSize,param.imgSize);
    ret.onlyatt=zeros(param.imgSize,param.imgSize);

    ret.otf=Airy_OtfWrite(ret.otf,ret,1,0,0);
    ret.onlyatt=getonlyatt(ret,0,0);
    ret.otfatt=ret.otf.*ret.onlyatt;
end