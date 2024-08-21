function k2fa = peakCCidx(obj,fAo0,fAp0,Zmask)
            
            Zo = double( obj.Ro<obj.Kotf );
            %tic
            %Tbegin = toc;
            Zn = fftshift( ifft2( fft2(Zo).*conj(fft2(Zo)) ) );
            CC = fftshift( ifft2( fft2(fAo0.*Zo).*conj(fft2(fAp0.*Zo)) ) );

%             rho = abs(CC)./abs(Zn);
            rho = abs(CC)./(abs(Zn)+0.001); 
            temp = rho .*Zmask;
% figure; pcolor(log(abs(rho)+1)); shading interp
% figure; pcolor(Zmask); shading interp  
% figure; pcolor(temp); shading interp  
            [~, v] = max(temp(:));
            
            pmax = mod(v-1,obj.w) + 1; % row
            qmax = (v-pmax)/obj.w + 1; % column
            px0 = qmax - (obj.wo+1);
            py0 = pmax - (obj.wo+1);
            k2fa = [py0 px0]; % nearest pixel approximation

end