function func=csr1d_drift(k,s)
% use R. Bosch's CER impedance, Eq.(4) of PRST-AB 11, 090702 (2008)

global s_ele egamma_vec;

egamma=interp1(s_ele,egamma_vec,s);

iEnd=1;
iSum=0;

if (auxr(s) < 1e10)
    %fprintf('CSR drift calculated within dipole...\n'); pause;
    func=0.0;
else
    [Lb,rho]=find_Lb_and_rho_for_csr1d_drift(s);              % unit: cm
    downstream_s=find_downstream_s(s);                        % unit: cm
    
    if (Lb(1)==0) % no dipoles upstream
        tmp=0.0;
    else
        
        [r,~]=size(downstream_s);
        
        % Case D
        % R. Bosch, Eq.(24), PRST-AB 10, 050701 (2007)
        for m=1:1:r
            if (downstream_s(m) > abs(rho(m))^(2/3)*(2*pi/k)^(1/3))
                tmp01(m)=min(downstream_s(m),egamma^2/k);
                tmp03(m)=2/tmp01(m);
            else
            	tmp03(m)=0.0;
            end
        end
        
        %{
        % R. Bosch, Eq.(4)
        for m=1:1:r
            tmp01=min(downstream_s(m),egamma^2/k);
            tmp02=rho(m)^(2/3)*(2*pi/k)^(1/3);
            tmp03(m)=2*log(tmp01/tmp02)/downstream_s(m);
            if (tmp03(m) < 0.0)
                tmp03(m)=0.0;
            end
        end
        %}
        
        % Case C
        tmp04=exp(-1i*k*Lb.^2./(6*rho.^2).*(Lb+3*downstream_s));
        tmp05=Lb+2*downstream_s;
        tmp06=-4*tmp04./tmp05;
        
        if (iSum==1); tmp07=sum(tmp06); end
        if (iEnd==1); tmp07=tmp06(end); end
        
        tmp=tmp03(end)+tmp07;
        %tmp=tmp03(end);
        
    end
    func=tmp;
end

