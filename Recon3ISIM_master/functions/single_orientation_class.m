function obj = single_orientation_class(OTFo,Kotf,Snoisy,PSFe)
            % constructor
            obj.OTFo = OTFo;
            obj.PSFe =PSFe;
            obj.Kotf = Kotf; 
            obj.w = size(Snoisy,1); % = OTFo size
            obj.wo = obj.w/2;
            
            x = linspace(0,obj.w-1,obj.w);
            y = linspace(0,obj.w-1,obj.w);
            [X,Y] = meshgrid(x,y);
            obj.Cv = (X-obj.wo) + 1i*(Y-obj.wo);
            obj.Ro = abs(obj.Cv);
            
            % G is a notch-filter (determined heuristically)
            % for suppressing out-of-focus signal of off-center components 
%             Rg = obj.Ro.*(512/obj.w);
            Rg = obj.Ro.*(1/4); % 1/4 determined the range of notch-filter 
            obj.G = 1 - exp(-0.05*Rg.^1.2);
            % edge tapering raw SIM images and computing its FFT
            fSnoisy = zeros(size(Snoisy));
            for i = 1:7
                temp = edgetaper(Snoisy(:,:,i),obj.PSFe);
                fSnoisy(:,:,i) =  fftshift(fft2(temp));                
            end
            obj.fSnoisy = fSnoisy;
            clear temp Snoisy fSnoisy
            
end