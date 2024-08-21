function [wFilter2] = WienerFilterW2_2D_new(param)
     
    wFilter2 = WienerFilterW2_2D( param );

end

function [wFilter]=WienerFilterW2_2D( sp )
    wFilter.sp=sp;
    wFilter.wDenom=updateCache2(sp);
end

function [ wDenom ] = updateCache2( sp )
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

function [wDenom]=addWienerDenominator_2D( wd,sp,d,b)
    siz=size(wd);
    w=siz(2);
    h=siz(1);
    dir=sp.Dir(d);
    cyclMicron=sp.cyclesPerMicron;
    cnt=[h/2+1,w/2+1];

    x=1:w;
    y=1:h;
    [x,y]=meshgrid(x,y);

    rad1=hypot(x-cnt(2)-(b-1)*dir.px,y-cnt(1)-(b-1)*dir.py)*cyclMicron;
    rad2=hypot(x-cnt(2)+(b-1)*dir.px,y-cnt(1)+(b-1)*dir.py)*cyclMicron;

    otfVal1=abs(getOtfVval(sp.OtfProvider,rad1)).^2;   
    otfVal2=abs(getOtfVval(sp.OtfProvider,rad2)).^2;

    if b==1
        if sp.OtfProvider.attStrength==0
            otfVal1=otfVal1/2;
            otfVal2=otfVal2/2;
        else
            otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.05,1.0*sp.OtfProvider.attFWHM)/2;
            otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.05,1.0*sp.OtfProvider.attFWHM)/2;
        end
    else
            otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.15,1.0*sp.OtfProvider.attFWHM)/1.0;
            otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.15,1.0*sp.OtfProvider.attFWHM)/1.0;
    end
    wd=wd+otfVal1+otfVal2;
    wDenom=wd;
end
function [ val ]= getOtfVval(ret,cycl)
    mask=cycl>ret.cutoff;
    cycl(mask)=0; 

    pos=cycl./ret.cyclesPerMicron;
    cpos=pos+1;
    lpos=floor(cpos);
    hpos=ceil(cpos);
    f=(cpos-lpos); 

    retl=ret.vals(lpos).*(1-f);
    reth=ret.vals(hpos).*f;
    val=retl+reth;
    val(mask)=0;
end
function [va]= valAttenuation(dist,str,fwhm)
    va=1-str*(exp(-power(dist,2)/(power(0.5*fwhm,2))));
end