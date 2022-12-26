function func=csr1d_nur_drift(k_wave,s)

global s_ele egamma_vec;

small=1e-6;
egamma=interp1(s_ele,egamma_vec,s);

iEnd=1;
iSum=0; % if set iSum=1, remember to modify for-loop on line 22 m=1:1:r

imethod=2; ratio=1.0;

if (auxr(s) < 1e10)
    func=0.0;
else
    [Lb,rho]=find_Lb_and_rho_for_csr1d_drift(s);              % unit: cm
    downstream_s=find_downstream_s(s);                        % unit: cm
    
    if (Lb(1)==0) % no dipoles upstream
        tmp=0.0;
    else
        [r,~]=size(downstream_s);
        for m=r:1:r
            L=downstream_s(m);
            theta_m=Lb(m)/rho(m);
            
            if (imethod==1)
            tmp01=f_kernel(L,theta_m,rho(m),egamma)*exp(-1i*k_wave*dz(L,theta_m,rho(m),egamma));
            tmp02=f_kernel(L,0.0,rho(m),egamma)*exp(-1i*k_wave*dz(L,0.0,rho(m),egamma));
            tmp03=@(z) f_kernel(L,z,rho(m),egamma);
            tmp04=@(z) exp(-1i*k_wave*dz(L,z,rho(m),egamma));
            tmp05=@(z) dz_dtheta(L,z,rho(m),egamma);
            tmp06=@(z) tmp03(z).*tmp04(z).*tmp05(z);
            tmp07=quadgk(tmp06,small,theta_m); warning off;
            tmp08(m)=-4/rho(m)*(tmp01+1i*k_wave*tmp07);
            end
            
            if (imethod==2)
            tmp03=@(z) df_kernel_dtheta(L,z,rho(m),egamma);
            tmp04=@(z) exp(-1i*k_wave*dz(L,z,rho(m),egamma));
            tmp05=@(z) tmp03(z).*tmp04(z);
            tmp06=quadgk(tmp05,small,ratio*theta_m); warning off;
            tmp08(m)=-4/rho(m)*(tmp06);
            end
            
        end
        
        % Case C
        tmp09=exp(-1i*k_wave*Lb.^2./(6*rho.^2).*(Lb+3*downstream_s));
        tmp10=Lb+2*downstream_s;
        tmp11=-4*tmp09./tmp10;
        
        if (iSum==1); tmp12=sum(tmp11); end
        if (iEnd==1); tmp12=tmp11(end); end
        
        tmp=tmp08(end)+tmp12;
        %tmp=tmp08(end);
    end
    func=tmp;
    
end
