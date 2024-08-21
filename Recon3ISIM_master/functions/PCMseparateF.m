function [fDo,fDp,fDm,fDp1,fDm1,fDp2,fDm2] = PCMseparateF(obj)
            fS1noisy = obj.fSnoisy(:,:,1);
            fS2noisy = obj.fSnoisy(:,:,2);
            fS3noisy = obj.fSnoisy(:,:,3);
            fS4noisy = obj.fSnoisy(:,:,4);
            fS5noisy = obj.fSnoisy(:,:,5);
            fS6noisy = obj.fSnoisy(:,:,6);
            fS7noisy = obj.fSnoisy(:,:,7);
            phase = obj.phaseA;
            %% Separating the three frequency components
            MF = 1.0;
            %% Transformation Matrix
            M = 0.5*[1 0.5*MF*exp(-1i*phase(1)) 0.5*MF*exp(+1i*phase(1)) 0.5*MF*exp(-1i*phase(2)) 0.5*MF*exp(+1i*phase(2)) 0.5*MF*exp(-1i*phase(3)) 0.5*MF*exp(+1i*phase(3));
                   1 0.5*MF*exp(-1i*(phase(1)+2*pi/7)) 0.5*MF*exp(+1i*(phase(1)+2*pi/7)) 0.5*MF*exp(-1i*(phase(2)+2*2*pi/7)) 0.5*MF*exp(+1i*(phase(2)+2*2*pi/7)) 0.5*MF*exp(-1i*(phase(3)+3*2*pi/7)) 0.5*MF*exp(+1i*(phase(3)+3*2*pi/7));
                   1 0.5*MF*exp(-1i*(phase(1)+4*pi/7)) 0.5*MF*exp(+1i*(phase(1)+4*pi/7)) 0.5*MF*exp(-1i*(phase(2)+2*4*pi/7)) 0.5*MF*exp(+1i*(phase(2)+2*4*pi/7)) 0.5*MF*exp(-1i*(phase(3)+3*4*pi/7)) 0.5*MF*exp(+1i*(phase(3)+3*4*pi/7));
                   1 0.5*MF*exp(-1i*(phase(1)+6*pi/7)) 0.5*MF*exp(+1i*(phase(1)+6*pi/7)) 0.5*MF*exp(-1i*(phase(2)+2*6*pi/7)) 0.5*MF*exp(+1i*(phase(2)+2*6*pi/7)) 0.5*MF*exp(-1i*(phase(3)+3*6*pi/7)) 0.5*MF*exp(+1i*(phase(3)+3*6*pi/7));
                   1 0.5*MF*exp(-1i*(phase(1)+8*pi/7)) 0.5*MF*exp(+1i*(phase(1)+8*pi/7)) 0.5*MF*exp(-1i*(phase(2)+2*8*pi/7)) 0.5*MF*exp(+1i*(phase(2)+2*8*pi/7)) 0.5*MF*exp(-1i*(phase(3)+3*8*pi/7)) 0.5*MF*exp(+1i*(phase(3)+3*8*pi/7));
                   1 0.5*MF*exp(-1i*(phase(1)+10*pi/7)) 0.5*MF*exp(+1i*(phase(1)+10*pi/7)) 0.5*MF*exp(-1i*(phase(2)+2*10*pi/7)) 0.5*MF*exp(+1i*(phase(2)+2*10*pi/7)) 0.5*MF*exp(-1i*(phase(3)+3*10*pi/7)) 0.5*MF*exp(+1i*(phase(3)+3*10*pi/7));
                   1 0.5*MF*exp(-1i*(phase(1)+12*pi/7)) 0.5*MF*exp(+1i*(phase(1)+12*pi/7)) 0.5*MF*exp(-1i*(phase(2)+2*12*pi/7)) 0.5*MF*exp(+1i*(phase(2)+2*12*pi/7)) 0.5*MF*exp(-1i*(phase(3)+3*12*pi/7)) 0.5*MF*exp(+1i*(phase(3)+3*12*pi/7))];
            %% Separting the components
            %===========================================================
           Minv = inv(M);
           fDo = Minv(1,1)*fS1noisy + Minv(1,2)*fS2noisy + Minv(1,3)*fS3noisy + Minv(1,4)*fS4noisy + Minv(1,5)*fS5noisy + Minv(1,6)*fS6noisy + Minv(1,7)*fS7noisy;
           fDp = Minv(2,1)*fS1noisy + Minv(2,2)*fS2noisy + Minv(2,3)*fS3noisy + Minv(2,4)*fS4noisy + Minv(2,5)*fS5noisy + Minv(2,6)*fS6noisy + Minv(2,7)*fS7noisy;
           fDm = Minv(3,1)*fS1noisy + Minv(3,2)*fS2noisy + Minv(3,3)*fS3noisy + Minv(3,4)*fS4noisy + Minv(3,5)*fS5noisy + Minv(3,6)*fS6noisy + Minv(3,7)*fS7noisy;
           fDp1 = Minv(4,1)*fS1noisy + Minv(4,2)*fS2noisy + Minv(4,3)*fS3noisy + Minv(4,4)*fS4noisy + Minv(4,5)*fS5noisy + Minv(4,6)*fS6noisy + Minv(4,7)*fS7noisy;
           fDm1 = Minv(5,1)*fS1noisy + Minv(5,2)*fS2noisy + Minv(5,3)*fS3noisy + Minv(5,4)*fS4noisy + Minv(5,5)*fS5noisy + Minv(5,6)*fS6noisy + Minv(5,7)*fS7noisy;
           fDp2 = Minv(6,1)*fS1noisy + Minv(6,2)*fS2noisy + Minv(6,3)*fS3noisy + Minv(6,4)*fS4noisy + Minv(6,5)*fS5noisy + Minv(6,6)*fS6noisy + Minv(6,7)*fS7noisy;
           fDm2 = Minv(7,1)*fS1noisy + Minv(7,2)*fS2noisy + Minv(7,3)*fS3noisy + Minv(7,4)*fS4noisy + Minv(7,5)*fS5noisy + Minv(7,6)*fS6noisy + Minv(7,7)*fS7noisy;
end
