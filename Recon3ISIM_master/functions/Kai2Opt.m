function CCo = Kai2Opt(obj,phaseShift,FcS1aT,FcS2aT,FcS3aT,FcS4aT,FcS5aT,FcS6aT)

            phaseShift2 = phaseShift(1);
            phaseShift3 = phaseShift(2);
            phaseShift4 = phaseShift(3);
            phaseShift5 = phaseShift(4);
            phaseShift6 = phaseShift(5);
            phaseShift7 = phaseShift(6);

            phaseShift1 = 0;
            phase2 = exp(-1i*phaseShift2) - exp(-1i*phaseShift1);
            phase3 = exp(-1i*phaseShift3) - exp(-1i*phaseShift1);
            phase4 = exp(-1i*phaseShift4) - exp(-1i*phaseShift1);
            phase5 = exp(-1i*phaseShift5) - exp(-1i*phaseShift1);
            phase6 = exp(-1i*phaseShift6) - exp(-1i*phaseShift1);
            phase7 = exp(-1i*phaseShift7) - exp(-1i*phaseShift1);
            
            phase2_2 = exp(-1i*2*phaseShift2) - exp(-1i*phaseShift1);
            phase3_2 = exp(-1i*2*phaseShift3) - exp(-1i*phaseShift1);
            phase4_2 = exp(-1i*2*phaseShift4) - exp(-1i*phaseShift1);
            phase5_2 = exp(-1i*2*phaseShift5) - exp(-1i*phaseShift1);
            phase6_2 = exp(-1i*2*phaseShift6) - exp(-1i*phaseShift1);
            phase7_2 = exp(-1i*2*phaseShift7) - exp(-1i*phaseShift1);
            
            phase2_3 = exp(-1i*3*phaseShift2) - exp(-1i*phaseShift1);
            phase3_3 = exp(-1i*3*phaseShift3) - exp(-1i*phaseShift1);
            phase4_3 = exp(-1i*3*phaseShift4) - exp(-1i*phaseShift1);
            phase5_3 = exp(-1i*3*phaseShift5) - exp(-1i*phaseShift1);
            phase6_3 = exp(-1i*3*phaseShift6) - exp(-1i*phaseShift1);
            phase7_3 = exp(-1i*3*phaseShift7) - exp(-1i*phaseShift1);
            
            % Transformation Matrix
            M = [phase2 conj(phase2) phase2_2 conj(phase2_2) phase2_3 conj(phase2_3);
                phase3 conj(phase3) phase3_2 conj(phase3_2) phase3_3 conj(phase3_3)
                phase4 conj(phase4) phase4_2 conj(phase4_2) phase4_3 conj(phase4_3)
                phase5 conj(phase5) phase5_2 conj(phase5_2) phase5_3 conj(phase5_3)
                phase6 conj(phase6) phase6_2 conj(phase6_2) phase6_3 conj(phase6_3)
                phase7 conj(phase7) phase7_2 conj(phase7_2) phase7_3 conj(phase7_3)];

            % Separating the components
            %===========================================================
            Minv = inv(M);
            FiSMap1 = Minv(1,1)*FcS1aT + Minv(1,2)*FcS2aT+Minv(1,3)*FcS3aT+Minv(1,4)*FcS4aT+Minv(1,5)*FcS5aT+Minv(1,6)*FcS6aT;
            FiSMam1 = Minv(2,1)*FcS1aT + Minv(2,2)*FcS2aT+Minv(2,3)*FcS3aT+Minv(2,4)*FcS4aT+Minv(2,5)*FcS5aT+Minv(2,6)*FcS6aT;
            FiSMap2 = Minv(3,1)*FcS1aT + Minv(3,2)*FcS2aT+Minv(3,3)*FcS3aT+Minv(3,4)*FcS4aT+Minv(3,5)*FcS5aT+Minv(3,6)*FcS6aT;
            FiSMam2 = Minv(4,1)*FcS1aT + Minv(4,2)*FcS2aT+Minv(4,3)*FcS3aT+Minv(4,4)*FcS4aT+Minv(4,5)*FcS5aT+Minv(4,6)*FcS6aT;
            FiSMap3 = Minv(5,1)*FcS1aT + Minv(5,2)*FcS2aT+Minv(5,3)*FcS3aT+Minv(5,4)*FcS4aT+Minv(5,5)*FcS5aT+Minv(5,6)*FcS6aT;
            FiSMam3 = Minv(6,1)*FcS1aT + Minv(6,2)*FcS2aT+Minv(6,3)*FcS3aT+Minv(6,4)*FcS4aT+Minv(6,5)*FcS5aT+Minv(6,6)*FcS6aT;
            %{
            figure;
            surf(log(abs(FiSMap)),'EdgeColor','none')
            colormap jet
            figure;
            surf(log(abs(FiSMam)),'EdgeColor','none')
            colormap jet
            kk
            %}

            CCo = abs( sum(sum( FiSMap1.*conj(FiSMam1).*obj.G )) )+...
                abs( sum(sum( FiSMap2.*conj(FiSMam2).*obj.G )) )+...
                abs( sum(sum( FiSMap3.*conj(FiSMam3).*obj.G )) );
end