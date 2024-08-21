function [FiSMao,FiSMap,FiSMam,FiSMap1,FiSMam1,FiSMap2,FiSMam2] = SeparatedComponents2D(...
             phaseShift,phaseShift0,FcS1aT,FcS2aT,FcS3aT,FcS4aT,FcS5aT,FcS6aT,FcS7aT)
% Aim: Unmixing the frequency components of raw SIM images
%   phaseShift,phaseShift0: illumination phase shifts
%   FcS1aT,FcS2aT,FcS3aT: FT of raw SIM images
%   FiSMao,FiSMap,FiSMam: unmixed frequency components of raw SIM images
               phaseShift1 = phaseShift(1);
               phaseShift2 = phaseShift(2);
               phaseShift3 = phaseShift(3);
               phaseShift4 = phaseShift(4);
               phaseShift5 = phaseShift(5);
               phaseShift6 = phaseShift(6);
               MF = 1.0;
%% Transformation Matrix
          M = 0.5*[1 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0) 0.5*MF*exp(-1i*phaseShift0) 0.5*MF*exp(+1i*phaseShift0);
                   1 0.5*MF*exp(-1i*phaseShift1) 0.5*MF*exp(+1i*phaseShift1) 0.5*MF*exp(-1i*2*phaseShift1) 0.5*MF*exp(+1i*2*phaseShift1) 0.5*MF*exp(-1i*3*phaseShift1) 0.5*MF*exp(+1i*3*phaseShift1);
                   1 0.5*MF*exp(-1i*phaseShift2) 0.5*MF*exp(+1i*phaseShift2) 0.5*MF*exp(-1i*2*phaseShift2) 0.5*MF*exp(+1i*2*phaseShift2) 0.5*MF*exp(-1i*3*phaseShift2) 0.5*MF*exp(+1i*3*phaseShift2);
                   1 0.5*MF*exp(-1i*phaseShift3) 0.5*MF*exp(+1i*phaseShift3) 0.5*MF*exp(-1i*2*phaseShift3) 0.5*MF*exp(+1i*2*phaseShift3) 0.5*MF*exp(-1i*3*phaseShift3) 0.5*MF*exp(+1i*3*phaseShift3);
                   1 0.5*MF*exp(-1i*phaseShift4) 0.5*MF*exp(+1i*phaseShift4) 0.5*MF*exp(-1i*2*phaseShift4) 0.5*MF*exp(+1i*2*phaseShift4) 0.5*MF*exp(-1i*3*phaseShift4) 0.5*MF*exp(+1i*3*phaseShift4);
                   1 0.5*MF*exp(-1i*phaseShift5) 0.5*MF*exp(+1i*phaseShift5) 0.5*MF*exp(-1i*2*phaseShift5) 0.5*MF*exp(+1i*2*phaseShift5) 0.5*MF*exp(-1i*3*phaseShift5) 0.5*MF*exp(+1i*3*phaseShift5);
                   1 0.5*MF*exp(-1i*phaseShift6) 0.5*MF*exp(+1i*phaseShift6) 0.5*MF*exp(-1i*2*phaseShift6) 0.5*MF*exp(+1i*2*phaseShift6) 0.5*MF*exp(-1i*3*phaseShift6) 0.5*MF*exp(+1i*3*phaseShift6)];

%% Separting the components
%===========================================================
          Minv = inv(M);

           FiSMao = Minv(1,1)*FcS1aT + Minv(1,2)*FcS2aT + Minv(1,3)*FcS3aT + Minv(1,4)*FcS4aT + Minv(1,5)*FcS5aT + Minv(1,6)*FcS6aT + Minv(1,7)*FcS7aT;
           FiSMap = Minv(2,1)*FcS1aT + Minv(2,2)*FcS2aT + Minv(2,3)*FcS3aT + Minv(2,4)*FcS4aT + Minv(2,5)*FcS5aT + Minv(2,6)*FcS6aT + Minv(2,7)*FcS7aT;
           FiSMam = Minv(3,1)*FcS1aT + Minv(3,2)*FcS2aT + Minv(3,3)*FcS3aT + Minv(3,4)*FcS4aT + Minv(3,5)*FcS5aT + Minv(3,6)*FcS6aT + Minv(3,7)*FcS7aT;
           FiSMap1 = Minv(4,1)*FcS1aT + Minv(4,2)*FcS2aT + Minv(4,3)*FcS3aT + Minv(4,4)*FcS4aT + Minv(4,5)*FcS5aT + Minv(4,6)*FcS6aT + Minv(4,7)*FcS7aT;
           FiSMam1 = Minv(5,1)*FcS1aT + Minv(5,2)*FcS2aT + Minv(5,3)*FcS3aT + Minv(5,4)*FcS4aT + Minv(5,5)*FcS5aT + Minv(5,6)*FcS6aT + Minv(5,7)*FcS7aT;
           FiSMap2 = Minv(6,1)*FcS1aT + Minv(6,2)*FcS2aT + Minv(6,3)*FcS3aT + Minv(6,4)*FcS4aT + Minv(6,5)*FcS5aT + Minv(6,6)*FcS6aT + Minv(6,7)*FcS7aT;
           FiSMam2 = Minv(7,1)*FcS1aT + Minv(7,2)*FcS2aT + Minv(7,3)*FcS3aT + Minv(7,4)*FcS4aT + Minv(7,5)*FcS5aT + Minv(7,6)*FcS6aT + Minv(7,7)*FcS7aT;

end 