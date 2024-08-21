function [OTFo,PSFe, Kotf, param] = Generate_OTF(Snoisy, param, NPixel)
%%  Generate approximate OTF
param.imgSize = size(Snoisy,1);
param.Size1 = NPixel; param.Size2 = NPixel;   
param.Iraw = Snoisy; 
param.micronsPerPixel = param.pixelsize;
param.cutoff=1000/(0.5*param.lambda/param.NA);                        

param.attStrength = 0;
param.cyclesPerMicron=1/(NPixel*param.micronsPerPixel);
param.OtfProvider=Airy_giveOtf(param,1);
param.sampleLateral=ceil(param.cutoff/param.cyclesPerMicron)+1;

OTFo = param.OtfProvider.otf;
OTFo = OTFo+eps; 
OTFo = OTFpost(OTFo);

PSFe = fspecial('gaussian',30,3.0); % for edgetapering

%% give the value of Kotf  
Kotf = OTFedgeF(OTFo); 
end
