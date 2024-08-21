function k2fa = TripeakCCidx(obj,fAo0,fAp0,fAp1,fAp2,Zmask)
            
            Zo = double( obj.Ro<obj.Kotf );
            %tic
            %Tbegin = toc;
            Zn = fftshift( ifft2( fft2(Zo).*conj(fft2(Zo)) ) );
            CC1 = fftshift( ifft2( fft2(fAo0.*Zo).*conj(fft2(fAp0.*Zo)) ) );
            CC2 = fftshift( ifft2( fft2(fAo0.*Zo).*conj(fft2(fAp1.*Zo)) ) );
            CC3 = fftshift( ifft2( fft2(fAo0.*Zo).*conj(fft2(fAp2.*Zo)) ) );
%             rho = abs(CC)./abs(Zn);
            rho1 = CC1./(abs(Zn)+0.001); 
            rho2 = CC2./(abs(Zn)+0.001); 
            rho3 = CC3./(abs(Zn)+0.001);
            rho = rho1+rho2+rho3;
            rho = abs(rho);
            temp = rho.*Zmask;
% figure; pcolor(log(abs(temp)+1)); shading interp
% figure; pcolor(Zmask); shading interp  
% figure; pcolor(temp); shading interp  
            [~, v] = max(temp(:));
            
            pmax = mod(v-1,obj.w) + 1; % row
            qmax = (v-pmax)/obj.w + 1; % column
            px0 = qmax - (obj.wo+1);
            py0 = pmax - (obj.wo+1);
            k2fa = [py0 px0]; % nearest pixel approximation

end
