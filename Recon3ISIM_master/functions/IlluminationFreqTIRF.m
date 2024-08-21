function k2fa1 = IlluminationFreqTIRF(obj,fDo,fDp,kA)
            % for noise suppression
            fAo0 = fDo.*obj.G.*(obj.Ro<obj.Kotf);
            fAp0 = fDp.*obj.G.*(obj.Ro<obj.Kotf);
            
            kv = kA(2) + 1i*kA(1);

            %Ro = abs(Cv);
            Rp = abs(obj.Cv-kv);
            Rm = abs(obj.Cv+kv);
            
            % Specify the search region for illumination freq peak
            Zmask = (1-(Rm<0.1*obj.Kotf ));%.*(1-(Rm<0.04*obj.Kotf )); % revised
%             k2fa = obj.peakCCidx(fAo0,fAp0,Zmask)
            k2fa = peakCCidx(obj,fAo0,fAp0,Zmask);
            
            % subpixel approximation            
            %----------------------------------
            if ( 2*obj.Kotf > obj.wo )
                t = 2*obj.w;
                fAo_temp = zeros(t,t);
                fAp_temp = zeros(t,t);
                fAo_temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = fAo0;
                fAp_temp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = fAp0;
                clear fAo fAp
                fAo0 = fAo_temp;
                fAp0 = fAp_temp;
                clear fAo_temp fAp_temp
            end
            %----------------------------------
            % the following parameters are defined here to reduce the
            % amount of computation within the Ifreq2optF function            
            t = size(fAo0,1);
            u = linspace(0,t-1,t);
            v = linspace(0,t-1,t);
            [U,V] = meshgrid(u,v);            
            Ap0 = ifft2(fAp0);
            %----------------------------------
            opt = 1;
%             Ifreq2opt0 = @(k2fa0)obj.Ifreq2optF(k2fa0,fAo0,Ap0,U,V,opt);
            Ifreq2opt0 = @(k2fa0)Ifreq2optF(obj,k2fa0,fAo0,Ap0,U,V,opt);
            options = optimset('LargeScale','off','Algorithm',...
                'active-set','MaxFunEvals',800,'MaxIter',800,'Display','notify');
            % options = optimset('LargeScale','off','Algorithm',...
            %             'active-set','MaxIter',200,'Display','iter');
            k2fa0 = k2fa;
            [k2fa1,fval] = fminsearch(Ifreq2opt0,k2fa0,options);
            %k2fa10
            %k2a = sqrt(k2fa1*k2fa1')
end