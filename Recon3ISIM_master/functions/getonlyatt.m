function [ onlyatt ] = getonlyatt( ret,kx,ky )
    w=ret.imgSize;
    h=ret.imgSize;
    siz=[h w];
    cnt=siz/2+1;
    kx=kx+cnt(2);
    ky=ky+cnt(1);
    onlyatt=zeros(h,w);

    y=1:h;
    x=1:w;
    [x,y]=meshgrid(x,y);
    rad=hypot(y-ky,x-kx);
    cycl=rad.*ret.cyclesPerMicron;
    onlyatt=valAttenuation(cycl,ret.attStrength,ret.attFWHM);

end