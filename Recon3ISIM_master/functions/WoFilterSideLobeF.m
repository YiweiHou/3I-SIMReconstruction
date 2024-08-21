function [FiSMaof,NoisePower] = WoFilterSideLobeF(obj,FiSMao,co,OBJsideP,SFo)

            OTFpower = obj.OTFo.*conj(obj.OTFo);

            % frequency beyond which NoisePower estimate to be computed
            NoiseFreq = obj.Kotf + 20;

            % NoisePower determination
            Zo = obj.Ro>NoiseFreq;
            nNoise = FiSMao.*Zo;
            NoisePower = sum(sum( nNoise.*conj(nNoise) ))./sum(sum(Zo));

            % Object Power determination
            OBJpower = OBJsideP.^2;

            %% Wiener Filtering
            FiSMaof = FiSMao.*(SFo.*conj(obj.OTFo)./NoisePower)./((SFo.^2).*OTFpower./NoisePower + co./OBJpower);
end 