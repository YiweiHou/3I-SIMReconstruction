function Esum = OBJparaOpt(obj,OBJpara0,LogF,coe)
            Ro = obj.Ro;
            Ro(obj.wo+1,obj.wo+1) = 1; % to avoid nan

            % approximated signal power calculation
            A = OBJpara0(1);
            B = OBJpara0(2);
            LogFe = A.*Ro + B;

            % range of frequency over which SSE is computed
%             Zloop = (Ro<0.75*obj.Kotf).*(Ro>0.12*obj.Kotf);
            if coe<=0.6
                Zloop = (Ro<1*obj.Kotf*coe).*(Ro>(1-coe)*obj.Kotf*coe);
            elseif coe>0.6 && coe<0.99
                Zloop = (Ro<1*obj.Kotf).*(Ro>0.3*obj.Kotf);
            else
                Zloop = (Ro<1*obj.Kotf).*(Ro>0.01*obj.Kotf);
            end
            % SSE computation
            Error = LogF - LogFe;
            Esum = sum(sum((Error.^2./Ro).*Zloop));
end
