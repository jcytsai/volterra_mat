function func=csr_shielding_pp(k,h,s)

% input: k(wave number) in cm^(-1), 
%        h(full height) in cm, and 
%        s(global longitudinal coordinate) in cm
% output: impedance (per unit length) Z in cgs unit

% Detailed formula can be found in Appendix of Agoh and
% Yokoya, PRSTAB 7 054403 (2004).

%{
c_speed=2.99792458e8;                      % unit: m/sec
Z0=376.73031;                              % unit: Ohm

k_tmp=k*1e2;                               % unit: m^(-1)
freq=k_tmp*c_speed/(2*pi);                 % unit: Hz
h_tmp=h*1e-2;                              % unit: m
rho_tmp=abs(auxr(s))/100;                  % unit: m

pathToScript=fullfile(pwd,'run_csrImpedance_in_matlab.sh');
height=sprintf('%10.10f',h_tmp);
radius=sprintf('%10.10f',rho_tmp);
Maxfreq=sprintf('%10.10f',10*freq);

cmdStr=[pathToScript ' ' height ' ' radius ' ' Maxfreq];
system(cmdStr);

load('Z_CSR_pp.txt');

%figure(1); plot(Z_CSR_pp(:,1),Z_CSR_pp(:,2),'r'); xlabel('k (m^{-1})'); ylabel('Re [Z] (Ohm)');
%figure(2); plot(Z_CSR_pp(:,1),Z_CSR_pp(:,3),'r'); xlabel('k (m^{-1})'); ylabel('Im [Z] (Ohm)');

k_ref=Z_CSR_pp(:,1)/1e2;                                 % unit: cm^(-1)
ZReal=Z_CSR_pp(:,2)*(4*pi/Z0)/(2*pi*rho_tmp*1e2);        % unit: rad
ZImag=-Z_CSR_pp(:,3)*(4*pi/Z0)/(2*pi*rho_tmp*1e2);       % unit: rad

ReZ=interp1(k_ref,ZReal,k);
ImZ=interp1(k_ref,ZImag,k);

func=(ReZ+1i*ImZ);
%}

p_max=200;
rho=abs(auxr(s));                      % unit: cm
tmp00=8*pi^2/h*(2/(k*rho))^(1/3);
tmp_sum=0.0;
for p=0:1:p_max
	beta_p=(2*p+1)*pi/h*(rho/(2*k^2))^(1/3);
	tmp_sum=tmp_sum+F0(beta_p);
end
func=tmp00*tmp_sum;

