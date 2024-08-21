function [va]=valIdealOTF(dist)
    if dist<0 || dist>1
        va=0;
        return;
    end
    va=(1/pi)*(2*acos(dist)-sin(2*acos(dist)));
end