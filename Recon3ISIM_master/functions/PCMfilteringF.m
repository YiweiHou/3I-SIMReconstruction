function [FDp,FDm,NpDp,NpDm] = PCMfilteringF(obj,co)
            FDP(:,:,1) = obj.fDp; FDP(:,:,2) = obj.fDp1; FDP(:,:,3) = obj.fDp2; 
            FDM(:,:,1) = obj.fDm; FDM(:,:,2) = obj.fDm1; FDM(:,:,3) = obj.fDm2; 
            for i = 1:3
            % suppressing out-of-focus signal of off-center components using
            fDp = FDP(:,:,i);
            fDm = FDM(:,:,i);
            % Object power parameters determination
            Aobj = obj.OBJpara(1);
            Bobj = obj.OBJpara(2);
            %% Duplex power (default)
            kv = obj.kAmean(i,2) + 1i*obj.kAmean(i,1); % vector along illumination direction
            Rp = abs(obj.Cv-kv);
            Rm = abs(obj.Cv+kv);
            OBJp = exp(Aobj.*Rp + Bobj);
            OBJm = exp(Aobj.*Rm + Bobj);

            k3 = round(obj.kAmean(i,:));
            OBJp(obj.wo+1+k3(1),obj.wo+1+k3(2)) = 0.25*OBJp(obj.wo+2+k3(1),obj.wo+1+k3(2))...
                + 0.25*OBJp(obj.wo+1+k3(1),obj.wo+2+k3(2))...
                + 0.25*OBJp(obj.wo+0+k3(1),obj.wo+1+k3(2))...
                + 0.25*OBJp(obj.wo+1+k3(1),obj.wo+0+k3(2));
            OBJm(obj.wo+1-k3(1),obj.wo+1-k3(2)) = 0.25*OBJm(obj.wo+2-k3(1),obj.wo+1-k3(2))...
                + 0.25*OBJm(obj.wo+1-k3(1),obj.wo+2-k3(2))...
                + 0.25*OBJm(obj.wo+0-k3(1),obj.wo+1-k3(2))...
                + 0.25*OBJm(obj.wo+1-k3(1),obj.wo+0-k3(2));

            % Filtering side lobes (off-center frequency components)
            SFo = obj.modFac(i);
            [fDpf,npDp] = WoFilterSideLobeF(obj,fDp,co,OBJm,SFo);
            [fDmf,npDm] = WoFilterSideLobeF(obj,fDm,co,OBJp,SFo);
            %% doubling Fourier domain size if necessary
%             if ( 2*obj.Kotf > obj.wo )
            if 1
                t = 2*obj.w;
                to = t/2;
                u = linspace(0,t-1,t);
                v = linspace(0,t-1,t);
                [U,V] = meshgrid(u,v);
                fDoTemp = zeros(t,t);
                fDpTemp = zeros(t,t);
                fDmTemp = zeros(t,t);
                fDoTemp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = obj.fDof;
                fDpTemp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = fDpf;
                fDmTemp(obj.wo+1:obj.w+obj.wo,obj.wo+1:obj.w+obj.wo) = fDmf;
                clear fDpf fDmf 
                fDof = fDoTemp;
                fDpf = fDpTemp;
                fDmf = fDmTemp;
                clear fDoTemp fDpTemp fDmTemp 
            else
                t = obj.w;
                to = t/2;
                u = linspace(0,t-1,t);
                v = linspace(0,t-1,t);
                [U,V] = meshgrid(u,v);
            end

            % Shifting the off-center frequency components to their correct location
            fDp1 = fft2(ifft2(fDpf).*exp( +1i.*2*pi*(obj.kAmean(i,2)/t.*(U-to) + obj.kAmean(i,1)/t.*(V-to)) ));
            fDm1 = fft2(ifft2(fDmf).*exp( -1i.*2*pi*(obj.kAmean(i,2)/t.*(U-to) + obj.kAmean(i,1)/t.*(V-to)) ));

           %% Shift induced phase error correction
            Cv = (U-to) + 1i*(V-to);
            Ro = abs(Cv);
            Rp = abs(Cv-kv);
            k2 = sqrt(obj.kAmean(i,:)*obj.kAmean(i,:)');

            % frequency range over which corrective phase is determined
            Zmask = (Ro < 0.8*k2).*(Rp < 0.8*k2);

            % corrective phase
            Angle0 = angle( sum(sum( fDof.*conj(fDp1).*Zmask )) );

            % phase correction
            fDp2 = exp(+1i*Angle0).*fDp1;
            fDm2 = exp(-1i*Angle0).*fDm1;
            FDp(:,:,i) = fDp2;
            FDm(:,:,i) = fDm2;
            NpDp(:,:,i) = npDp;
            NpDm(:,:,i) = npDm;
            end
end