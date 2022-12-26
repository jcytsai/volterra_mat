function Zk_vs_s(lambda_start01,lambda_end01,scan_num01,s_array)

% lambda_start01, lambda_start01: unit, cm
% s_array, s_ele: unit: cm
% NOTE: bending radius will be calculated by function auxr(s)

global s_ele C_ele sig_x_ele sig_y_ele egamma_vec
global iCSR_ss iCSR_tr iCSR_drift iLSC

factor=120*pi/(4*pi);

if (scan_num01==1)
    lambda=lambda_start01;
    for m=1:1:length(s_array)
        
        %if (m==1492)
        %    pause;
        %end
        
        s=s_array(m);
        Cpf=interp1(s_ele,C_ele,s);
        egamma=interp1(s_ele,egamma_vec,s);
        Sx=interp1(s_ele,sig_x_ele,s);
        Sy=interp1(s_ele,sig_y_ele,s);
        %rb=0.5*1.7*(Sx+Sy);
        lambda_s(m)=lambda/Cpf;
        k_s=2*pi/lambda_s(m);
        
        if (iCSR_ss==1)
            Z_CSR_ss(m)=csr1d_nur(k_s,s,egamma,iCSR_ss);   % NUR model
            %Z_CSR_ss(m)=csr1d(k_s,s,iCSR_ss);             % UR model
        else
            Z_CSR_ss(m)=0.0;
        end
        
        if (iCSR_tr==1)
            Z_CSR_tr(m)=csr1d_tr(k_s,s);
        else
            Z_CSR_tr(m)=0.0;
        end
        
        if (iCSR_drift==1)
            %Z_CSR_drif(m)=csr1d_drift(k_s,s);         % Bosch model
            Z_CSR_drif(m)=csr1d_nur_drift(k_s,s);      % Li model
        else
            Z_CSR_drif(m)=0.0;
        end
        
        if (iLSC==1)
            Z_LSC(m)=lsc1d(k_s,Sx,Sy,s,1);             % on-axis model
        else
            Z_LSC(m)=0.0;
        end
    end
    Z_CSR=Z_CSR_ss+Z_CSR_tr+Z_CSR_drif;
    Z_tot=Z_CSR+Z_LSC;
    
    figure(201); plot(s_array/100,factor*abs(Z_tot),'k'); xlabel('s (m)'); ylabel('Z(k;s) (\Omega/cm)'); title('CSR+LSC impedance'); hold on;
    figure(201); plot(s_array/100,factor*real(Z_tot),'r'); hold on;
    figure(201); plot(s_array/100,factor*imag(Z_tot),'b'); hold on;
    figure(202); plot(s_array/100,factor*abs(Z_CSR),'k'); xlabel('s (m)'); ylabel('|Z(k;s)| (\Omega/cm)'); title('CSR impedance'); hold on;
    figure(202); plot(s_array/100,factor*real(Z_CSR),'r'); hold on;
    figure(202); plot(s_array/100,factor*imag(Z_CSR),'b'); hold on;
    figure(203); plot(s_array/100,factor*abs(Z_LSC),'k'); xlabel('s (m)'); ylabel('|Z(k;s)| (\Omega/cm)'); title('LSC impedance'); hold on;
    figure(204); plot(s_array/100,lambda_s*1e4,'k'); xlabel('s (m)'); ylabel('\lambda(s) (\mum)');
else
    lambda_array=linspace(lambda_start01,lambda_end01,scan_num01);
    for m=1:1:length(s_array)
        s=s_array(m);
        Cpf=interp1(s_ele,C_ele,s);
        egamma=interp1(s_ele,egamma_vec,s);
        Sx=interp1(s_ele,sig_x_ele,s);
        Sy=interp1(s_ele,sig_y_ele,s);
        
        for n=1:1:scan_num01
        lambda_s=lambda_array(n)/Cpf;
        k_s=2*pi/lambda_s;
        
        if (iCSR_ss==1)
            Z_CSR_ss(m,n)=csr1d_nur(k_s,s,egamma,iCSR_ss);   % NUR model
            %Z_CSR_ss(m)=csr1d(k_s,s,iCSR_ss);             % UR model
        else
            Z_CSR_ss(m,n)=0.0;
        end
        
        if (iCSR_tr==1)
            Z_CSR_tr(m,n)=csr1d_tr(k_s,s);
        else
            Z_CSR_tr(m,n)=0.0;
        end
        
        if (iCSR_drift==1)
            %Z_CSR_drif(m,n)=csr1d_drift(k_s,s);         % Bosch model
            Z_CSR_drif(m,n)=csr1d_nur_drift(k_s,s);      % Li model
        else
            Z_CSR_drif(m,n)=0.0;
        end
        
        if (iLSC==1)
            Z_LSC(m,n)=lsc1d(k_s,Sx,Sy,s,1);             % on-axis model
        else
            Z_LSC(m,n)=0.0;
        end
        end
    end
    Z_CSR=Z_CSR_ss+Z_CSR_tr+Z_CSR_drif;
    Z_tot=Z_CSR+Z_LSC;
    
    figure(201); mesh(lambda_array*1e4,s_array/100,factor*abs(Z_tot)); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('Z(k;s) (\Omega/cm)'); title('CSR+LSC impedance'); hold on;
    %figure(201); mesh(lambda_array*1e4,s_array/100,factor*real(Z_tot)); hold on;
    %figure(201); mesh(lambda_array*1e4,s_array/100,factor*imag(Z_tot)); hold on;
    figure(202); mesh(lambda_array*1e4,s_array/100,factor*abs(Z_CSR)); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('|Z(k;s)| (\Omega/cm)'); title('CSR impedance'); hold on;
    %figure(202); mesh(lambda_array*1e4,s_array/100,factor*real(Z_CSR)); hold on;
    %figure(202); mesh(lambda_array*1e4,s_array/100,factor*imag(Z_CSR)); hold on;
    figure(203); mesh(lambda_array*1e4,s_array/100,factor*abs(Z_LSC)); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('|Z(k;s)| (\Omega/cm)'); title('LSC impedance'); hold on;
end