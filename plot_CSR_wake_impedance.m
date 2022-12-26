%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is independent of the main program {volterra_mat_for_GUI.m}
% and intends to illustrate various free-space CSR impedance models, 
% including steady-state, entrance transient and exit transients (drift CSR
% and "Case C"). It can also calculate the wakefield for a given arbitrary
% bunch distribution. (to be implemented)
% 
% This program uses cgs units.
%
% Program written by Cheng-Ying Tsai, jcytsai@vt.edu
% v1: March 21, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
format long

const_A=3^(-1/3)*gamma(2/3)*(1i*sqrt(3)-1);

iss=1;
itr=0;
idrif=1;


rho=2.0e2;                                   % bending radius in cm
Lb=20;                                       % dipole length in cm
egamma=0.75e3/0.511;                         % beam energy in GeV


lambda_start=0.1e-4;                           % unit: cm
lambda_end=100e-4;                           % unit: cm
scan_num=2000;
lambda_array=linspace(lambda_start,lambda_end,scan_num);
k_array=2*pi./lambda_array;

% steady-state CSR impedance
if (iss==1)
    Z_CSR_ss=-1i*k_array.^(1/3)*const_A*rho^(-2/3);
else
    Z_CSR_ss=0.0;
end


% entrance transient CSR impedance
if (itr==1)
    local_s=10;                               % local coordinate measured from dipole entrance, in cm
    z_L=local_s^3/(24*rho^2);

    % D. Zhou's version
    mu=k_array*z_L;
    tmp01=-4/local_s.*exp(-1i*4*mu);
    tmp02=4/(3*local_s)*(1i*mu).^(1/3).*mfun('GAMMA',-1/3,1i*mu);
    Z_CSR_tr=tmp01+tmp02;
else
    Z_CSR_tr=0.0;
end



% exit transient CSR impedances
if (idrif==1)
    downstream_s=100;                         % local coordinate measured from dipole exit, in cm        
    
    % R. Bosch, Eq.(24), PRST-AB 10, 050701 (2007)
    for m=1:1:length(k_array)
        k=k_array(m);
        if (downstream_s > rho^(2/3)*(2*pi/k)^(1/3))
        	tmp01(m)=min(downstream_s,egamma^2/k);
            tmp03(m)=2/tmp01(m);
        else
            tmp03(m)=0.0;
        end
        
        tmp04=exp(-1i*k*Lb^2/(6*rho^2)*(Lb+3*downstream_s));
        tmp05=Lb+2*downstream_s;
        tmp06(m)=-4*tmp04/tmp05;
        
        tmp(m)=tmp03(m)+tmp06(m);
        
    end
    Z_CSR_drif=tmp;
else
    Z_CSR_drif=0.0;
end

Z_tot=Z_CSR_ss+Z_CSR_tr+Z_CSR_drif;

% plots
if (iss==1); figure(1); plot(lambda_array*1e4,real(Z_CSR_ss),'rs-'); xlabel('\lambda (\mum)'); ylabel('Real part of impedances'); grid off; hold on; end
if (itr==1); figure(1); plot(lambda_array*1e4,real(Z_CSR_tr),'gd-'); xlabel('\lambda (\mum)'); ylabel('Real part of impedances'); grid off; hold on; end
if (idrif==1); figure(1); plot(lambda_array*1e4,real(Z_CSR_drif),'b-'); xlabel('\lambda (\mum)'); ylabel('Real part of impedances'); grid off; hold on; end

if (iss==1); figure(2); plot(lambda_array*1e4,imag(Z_CSR_ss),'rs-'); xlabel('\lambda (\mum)'); ylabel('Imaginary part of impedances'); grid off; hold on; end
if (itr==1); figure(2); plot(lambda_array*1e4,imag(Z_CSR_tr),'gd-'); xlabel('\lambda (\mum)'); ylabel('Imaginary part of impedances'); grid off; hold on; end
if (idrif==1); figure(2); plot(lambda_array*1e4,imag(Z_CSR_drif),'b-'); xlabel('\lambda (\mum)'); ylabel('Imaginary part of impedances'); grid off; hold on; end

figure(3); plot(lambda_array*1e4,real(Z_tot),'r'); xlabel('\lambda (\mum)'); ylabel('Real part of total impedances'); grid off; hold on;
figure(3); plot(lambda_array*1e4,real(Z_CSR_ss),'r--'); xlabel('\lambda (\mum)'); ylabel('Real part of total impedances'); grid off; hold on;

figure(4); plot(lambda_array*1e4,imag(Z_tot),'b'); xlabel('\lambda (\mum)'); ylabel('Imaginary part of total impedances'); grid off; hold on;
figure(4); plot(lambda_array*1e4,imag(Z_CSR_ss),'b--'); xlabel('\lambda (\mum)'); ylabel('Imaginary part of total impedances'); grid off; hold on;



