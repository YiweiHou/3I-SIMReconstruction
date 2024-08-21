function [fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2] = PCMseparate_tempF(obj)
            
            fS1aTnoisy = obj.fSnoisy(:,:,1);
            fS2aTnoisy = obj.fSnoisy(:,:,2);
            fS3aTnoisy = obj.fSnoisy(:,:,3);
            fS4aTnoisy = obj.fSnoisy(:,:,4);
            fS5aTnoisy = obj.fSnoisy(:,:,5);
            fS6aTnoisy = obj.fSnoisy(:,:,6);
            fS7aTnoisy = obj.fSnoisy(:,:,7);
            
            fDuplex2 = fS2aTnoisy - fS1aTnoisy;
            fDuplex3 = fS3aTnoisy - fS1aTnoisy;
            fDuplex4 = fS4aTnoisy - fS1aTnoisy;
            fDuplex5 = fS5aTnoisy - fS1aTnoisy;
            fDuplex6 = fS6aTnoisy - fS1aTnoisy;
            fDuplex7 = fS7aTnoisy - fS1aTnoisy;
            
            %% Optimizing the relative phases (computational time ~1.9s)
%             Kai2Opt0 = @(phase0)obj.Kai2Opt(phase0,fDuplex2,fDuplex3);
            Kai2Opt0 = @(phase0)Kai2Opt(obj,phase0,fDuplex2,fDuplex3,fDuplex4,fDuplex5,fDuplex6,fDuplex7);
            
            options = optimset('LargeScale','off','Algorithm',...
                'active-set','MaxFunEvals',800,'MaxIter',800,'Display','notify');
            phase0 = [2*pi/7 4*pi/7 6*pi/7 8*pi/7 10*pi/7 12*pi/7]; % initial guess
           % [phaseShift,fval] = fminsearch(Kai2Opt0,phase0,options);
            % Separating the three frequency components
            phaseShift = phase0;
 phaseShift0=0;
            [fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2] = SeparatedComponents2D(...
                phaseShift,phaseShift0,fS1aTnoisy,fS2aTnoisy,fS3aTnoisy,fS4aTnoisy,fS5aTnoisy,fS6aTnoisy,fS7aTnoisy);   
end