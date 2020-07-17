function peakTangent=givepeaktangent(obj)
        
        %Non-dimension shear and normal components
        dshear=obj.Disp(1);
        lcr=obj.Lambdacr;
        dcr=obj.CriticalDisp;
        scr=obj.PeakTraction;
        dnormal=lcr*dcr*1.01;           % dnormal right passes the lcr*dcr
        nds=dshear/dcr;
        ndn=dnormal/dcr;
        lambdae=sqrt(nds^2+ndn^2);      % nondimensional effective displacement
        % Assume zero crack stiffness if crack opening is zero, however, a
        % sudden jump in crack stiffness will be observed when lamdae
        % increases from zero to a tiny value. NEED TO BE DISCUSSED LATER
        % 02/04/2019. Maybe it is the reason why Khoei in 2010 used
        % bilinear cohesive law is to avoid such discontinuity.
        if lambdae>=obj.LambdaeL
            % Loading
            obj.LambdaeL=lambdae;
            obj.Loading=1;
            temp1=dcr*scr/(1-lcr);
            Tss=-temp1*(dshear/lambdae/dcr^2)^2+...
                (1-lambdae)*temp1*(1/lambdae/dcr^2-dshear^2/lambdae^3/dcr^4);
            Tsn=-temp1/lambdae^3*dshear*dnormal/dcr^4;
            Tnn=-temp1*(dnormal/lambdae/dcr^2)^2+...
                (1-lambdae)*temp1*(1/lambdae/dcr^2-dnormal^2/lambdae^3/dcr^4);
        else
            % unloading
            obj.Loading=-1;
            Tss=scr/dcr*(1-obj.LambdaeL)/(1-lcr)/obj.LambdaeL;
            Tsn=0;
            Tnn=Tss;
        end
        peakTangent=[Tss,Tsn;Tsn,Tnn];
    end