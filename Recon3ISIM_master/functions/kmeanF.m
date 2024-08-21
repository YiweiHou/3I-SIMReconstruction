function [kmean] =  kmeanF(obj)  
            kA = zeros(3,2);
            for i = 1:3
                fSa1Tnoisy0 = obj.sepf(:,:,1);
                fSa1Tnoisyp = obj.sepf(:,:,i+1);
                kA(i,:) = IlluminationFreqTIRF(obj,fSa1Tnoisy0,fSa1Tnoisyp,obj.kAo(i,:)); 
            end
            kmean = kA;
            %kmag = sqrt(kmean*kmean');
end