function [fDof,NoisePower] = WoFilterCenterF(obj,co)            

            % OTF Power 
            OTFpower = obj.OTFo.*conj(obj.OTFo);

            % frequency beyond which NoisePower estimate to be computed
            NoiseFreq = obj.Kotf + 20;

            % NoisePower determination
            Zo = obj.Ro>NoiseFreq;
            nNoise = obj.fDo.*Zo;
            NoisePower = sum(sum( nNoise.*conj(nNoise) ))./sum(sum(Zo));

            % Object Power determination
            A = obj.OBJpara(1);
            B = obj.OBJpara(2);
            LogFe = A*obj.Ro + B ;
            OBJo = exp(LogFe);
            OBJpower = OBJo.^2 - 1.0*NoisePower;

            %% Wiener Filtering
            SFo = 1;
            fDof = obj.fDo.*(SFo.*conj(obj.OTFo)./NoisePower)./((SFo.^2).*OTFpower./NoisePower + co./OBJpower);
end