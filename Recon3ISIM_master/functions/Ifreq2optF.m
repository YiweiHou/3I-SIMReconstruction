function [CCop] = Ifreq2optF(~,k2fa,fAoT,Ap0,U,V,opt)
            t = size(fAoT,1);
            to = t/2;
            
            ApT = exp( +1i.*2*pi*( k2fa(2)/t.*(U-to)+k2fa(1)/t.*(V-to) ) ).*Ap0;
            fApT0 = fft2( ApT );
            
            mA = sum(sum( fAoT.*conj(fApT0) ));
            mA = mA./sum(sum( fApT0.*conj(fApT0) ));

            if opt>0
                CCop = -abs(mA);
            else
                CCop = angle(mA);
            end
            
end