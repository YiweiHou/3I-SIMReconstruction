function [obj] = KapproxEstimationF(obj,flag)
%             [fDo_temp,fDp_temp,~] = obj.PCMseparate_tempF();
            [fDo_temp,fDp_temp,~,fDp_temp1,~,fDp_temp2,~] = PCMseparate_tempF(obj);
%  figure; pcolor(log(abs(fDo_temp)+1)); shading interp  
%  figure; pcolor(log(abs(fDp_temp)+1)); shading interp 
%  figure; pcolor(log(abs(fDp_temp1)+1)); shading interp 
%  figure; pcolor(log(abs(fDp_temp2)+1)); shading interp 
 
            % for noise suppression
            fAo0 = fDo_temp.*obj.G.*(obj.Ro<obj.Kotf);
            fAp0 = fDp_temp.*obj.G.*(obj.Ro<obj.Kotf);            
            fAp1 = fDp_temp1.*obj.G.*(obj.Ro<obj.Kotf);
            fAp2 = fDp_temp2.*obj.G.*(obj.Ro<obj.Kotf); 
            % Specify the search region for illumination freq peak
            if ~flag 
                Zmask = 1.*(obj.Ro>0.92*obj.Kotf ).*(obj.Ro<1.5*obj.Kotf ); % Turf pattern Spread spectrum coef >1.9
            else
                Zmask = 1.*(obj.Ro>0.5*obj.Kotf ).*(obj.Ro<0.95*obj.Kotf );  % Spread spectrum coef <1.9
            end

            kAo = zeros(3,2);
            k2fa = TripeakCCidx(obj,fAo0,fAp0,fAp1,fAp2,Zmask);
            kAo(1,:) = peakCCidx(obj,fAo0,fAp0,Zmask);
            kAo(2,:) = peakCCidx(obj,fAo0,fAp1,Zmask);
            kAo(3,:) = peakCCidx(obj,fAo0,fAp2,Zmask); 
            obj.sepf(:,:,1) = fDo_temp;
            obj.sepf(:,:,2) = fDp_temp;
            obj.sepf(:,:,3) = fDp_temp1;
            obj.sepf(:,:,4) = fDp_temp2;    % 提取分离频谱，估计频率k值 
            obj.kAo = kAo;
end