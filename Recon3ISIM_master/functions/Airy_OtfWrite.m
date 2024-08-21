function [ ret ] = Airy_OtfWrite( vec,otf,band, kx, ky )
    ret=otfToVector(vec,otf,band,kx,ky,0,1);
end