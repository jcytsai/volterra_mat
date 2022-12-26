function func=csr1d_nur_drift(k_wave,s)

global s_ele egamma_vec;

small=1e-6;
egamma=interp1(s_ele,egamma_vec,s);

iEnd=1;
iSum=0; % if set iSum=1, remember to modify for-loop on line 22 m=1:1:r


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
            tmp01=f_kernel(L,theta_m,rho(m),egamma)*exp(-1i*k_wave*dz(L,theta_m,rho(m),egamma));
            tmp02=f_kernel(L,0.0,rho(m),egamma)*exp(-1i*k_wave*dz(L,0.0,rho(m),egamma));
            tmp03=@(z) f_kernel(L,z,rho(m),egamma);
            tmp04=@(z) exp(-1i*k_wave*dz(L,z,rho(m),egamma));
            tmp05=@(z) dz_dtheta(L,z,rho(m),egamma);
            tmp06=@(z) tmp03(z).*tmp04(z).*tmp05(z);
            tmp07=quadgk(tmp06,small,theta_m); warning off;
            tmp08(m)=4/rho(m)*(tmp01+1i*k_wave*tmp07);
        end
        tmp=tmp08(end);
    end
    func=tmp;
end
