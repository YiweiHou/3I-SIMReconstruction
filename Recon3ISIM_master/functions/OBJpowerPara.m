function OBJparaA = OBJpowerPara(obj,OTFo,coe)

            LogF = log(abs(obj.fDo)) - log(OTFo);

            %% object power parameters through optimization
             OBJparaOpt0 = @(OBJpara0)OBJparaOpt(obj,OBJpara0,LogF,coe);
            
            options = optimset('LargeScale','off','Algorithm',...
                'active-set','MaxFunEvals',800,'MaxIter',800,'Display','notify');

            % obtaining crude initial guesses for Aobj and Bobj 
            Zm = (obj.Ro>0.4*obj.Kotf).*(obj.Ro<0.5*obj.Kotf);
            Bobj = sum(sum(LogF.*Zm))./sum(sum(Zm));
            Aobj = -0.5;
            OBJpara0 = [Aobj Bobj];

            % optimization step
            [OBJparaA,fval] = fminsearch(OBJparaOpt0,OBJpara0,options);

end