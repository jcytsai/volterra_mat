% THIS PROGRAM CALCULATES SASE FEL AMPLIFIER. THE INITIAL SHOT-NOISE POWER
% IS ESTIMATED FOLLOWING MCNEIL.

%------------------ Read input values from GUI ---------------------------%

HH = findobj(gcf,'Tag','I_b_fin');          I_b_fin = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','sigma_delta_fin');  sigma_delta_fin = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','sig_x');            sig_x = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','emit_nx_fin');      emit_nx_fin = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','bk_ini');           bk_ini = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','Pk_ini');           Pk_ini = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','und_beta_ave');     und_beta_ave = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','und_lambda');       und_lambda = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','und_L');            und_L = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','lambda_r');         lambda_r = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','P_L');              P_L = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','MX');               MX = HH.Value;

%--------------------------- GUI-INPUT PARAMETERS ------------------------%
%{
MX=0;
und_lambda=0.035;                                   % unit: m
und_L=6*68*und_lambda;                              % unit: m
lambda_r=7.3e-9;                                    % unit: m
P_L=10;                                             % seed power, unit: MW
bk_ini=3e-4;                                        % ini. shot-noise dmod amp.
Pk_ini=1e-5;                                        % ini. shot-noise emod amp.
%}
%-------------------------------------------------------------------------%


%--------------------------- DERIVED PARAMETERS --------------------------%
%{
I_b_fin=I_b*C_ele(end);                             % unit: Amp
sigma_delta_fin=SES_fin/(egamma_vec(end)*mc2);
sig_x=sig_x_ele(end)/100;                           % unit: m
emit_nx_fin=emitx/100*egamma;                       % unit: m
if (iIBS==1)
    emit_nx_fin=emit_gx_IBS/100*egamma;
end
und_beta_ave=sig_x^2/(emit_nx_fin/egamma_vec(end));
%}
und_ku=2*pi/und_lambda;
und_kL=2*pi/lambda_r;
und_K=sqrt(2*(2*(egamma_vec(end))^2*lambda_r/und_lambda-1));
und_JJ=besselj(0,und_K^2/(4+2*und_K^2))-besselj(1,und_K^2/(4+2*und_K^2));
rho_und=2.1*(I_b_fin*und_K^2*und_JJ^2*(und_lambda*1e2)^2/((egamma_vec(end))^3*(sig_x*1e6)^2))^(1/3);
Lg_und=und_lambda/(4*pi*rho_und);

if (MX==1)
    Lg_und_1D=Lg_und;
    MX_factor=Lambda_MX(Lg_und,und_kL,und_lambda,sig_x,und_beta_ave,emit_nx_fin/egamma_vec(end),sigma_delta_fin);
    Lg_und=Lg_und_1D*(1+MX_factor);
end

A_hat_array=Gf_31*Pk_ini/rho_und;                    % scaled emod amp.
B_hat_array=Gf_11*bk_ini;                            % absolute dmod amp.

tau_und_end=und_L/Lg_und;
%{
tmp01=und_lambda*1e2*und_K*und_JJ;
tmp02=(egamma_vec(end))^2*rho_und^2;
tmp03=sqrt(P_L)/(sig_x*1e6);
E_0=1.095e4*tmp03;                                  % unit: MV/m
a0_und=6.03*tmp01/tmp02*tmp03;                      % normalized seed field
%}
N_lambda=I_b_fin*lambda_r/(c_speed/100*charge);     % number of electrons per rad. wavelength
tmp04=6*sqrt(pi)*rho_und;
tmp05=N_lambda*sqrt(log(N_lambda/rho_und));
a0_und=sqrt(tmp04/tmp05); 

b_H=0.0;                                            % dmod @ seed wavelength, assumed
P_H=0.0;                                            % emod @ seed wavelength, assumed

und_L_array=linspace(0,und_L,201);
%und_L_array=linspace(0,9*Lg_und,201);
PS_phi_array=linspace(0,0,201);
s_hat_array=und_L_array/Lg_und;

kapa_1=exp(1i*pi/6);
kapa_2=exp(1i*5*pi/6);
kapa_3=exp(-1i*pi/2);

C11=(kapa_2*kapa_3)/((kapa_1-kapa_2)*(kapa_1-kapa_3));
C12=(kapa_2+kapa_3)/((kapa_1-kapa_2)*(kapa_1-kapa_3));
C13=1i/rho_und/((kapa_1-kapa_2)*(kapa_1-kapa_3));
C21=(kapa_1*kapa_3)/((kapa_1-kapa_2)*(-kapa_2+kapa_3));
C22=(kapa_1+kapa_3)/((kapa_1-kapa_2)*(-kapa_2+kapa_3));
C23=1i/rho_und/((kapa_1-kapa_2)*(-kapa_2+kapa_3));
C31=(kapa_1*kapa_2)/((kapa_1-kapa_3)*(kapa_2-kapa_3));
C32=(kapa_1+kapa_2)/((kapa_1-kapa_3)*(kapa_2-kapa_3));
C33=1i/rho_und/((kapa_1-kapa_3)*(kapa_2-kapa_3));

% SET INITIAL STATE VECTOR
a_und=[a0_und zeros(1,length(PS_phi_array))];
b_und=[b_H zeros(1,length(PS_phi_array))];
P_und=[P_H zeros(1,length(PS_phi_array))];

ind_sat=length(b_und);
z_sat=s_hat_array(end);

for m=1:1:length(PS_phi_array)
    if (m==1)
        s=s_hat_array(m);
    else
        s=s_hat_array(m)-s_hat_array(m-1);
    end
    phi=PS_phi_array(m);
    M11=C11*exp(kapa_1*s)+C21*exp(kapa_2*s)+C31*exp(kapa_3*s);
    M12=C12*exp(kapa_1*s-1i*phi)+C22*exp(kapa_2*s-1i*phi)+C32*exp(kapa_3*s-1i*phi);
    M13=C13*exp(kapa_1*s+1i*phi)+C23*exp(kapa_2*s+1i*phi)+C33*exp(kapa_3*s+1i*phi);
    M21=-(kapa_1*C11*exp(kapa_1*s)+kapa_2*C21*exp(kapa_2*s)+kapa_3*C31*exp(kapa_3*s));
    M22=-(kapa_1*C12*exp(kapa_1*s-1i*phi)+kapa_2*C22*exp(kapa_2*s-1i*phi)+kapa_3*C32*exp(kapa_3*s-1i*phi));
    M23=-(kapa_1*C13*exp(kapa_1*s+1i*phi)+kapa_2*C23*exp(kapa_2*s+1i*phi)+kapa_3*C33*exp(kapa_3*s+1i*phi));
    M31=-1i*rho_und*(kapa_1^2*C11*exp(kapa_1*s+1i*phi)+kapa_2^2*C21*exp(kapa_2*s+1i*phi)+kapa_3^2*C31*exp(kapa_3*s+1i*phi));
    M32=-1i*rho_und*(kapa_1^2*C12*exp(kapa_1*s)+kapa_2^2*C22*exp(kapa_2*s)+kapa_3^2*C32*exp(kapa_3*s));
    M33=-1i*rho_und*(kapa_1^2*C13*exp(kapa_1*s+2*1i*phi)+kapa_2^2*C23*exp(kapa_2*s+2*1i*phi)+kapa_3^2*C33*exp(kapa_3*s+2*1i*phi));
    
    a_und(m+1)=M11*a_und(m)+M12*b_und(m)+M13*P_und(m);
    b_und(m+1)=M21*a_und(m)+M22*b_und(m)+M23*P_und(m);
    P_und(m+1)=M31*a_und(m)+M32*b_und(m)+M33*P_und(m);
    
    % DETERMINE SATURATION LOCATION
    if (abs(abs(b_und(m))-1)<0.1)
        ind_sat=m;
        %z_sat=und_L_array(ind_sat);
        z_sat=s_hat_array(ind_sat);
    end
end

tmp13=abs(a_und)*sig_x*1e6*(egamma_vec(end))^2*rho_und^2;
tmp14=6.03*und_lambda*1e2*und_K*und_JJ;
Power_und=(tmp13/tmp14).^2;                         % unit: MW
P_sat=Power_und(ind_sat);                           % unit: MW
%
if (length(PS_phi_array)>5)
    %figure(102); subplot(3,1,1); plot(und_L_array,abs(a_und(1:end-1)),'rs'); xlabel('Radiator s (m)'); ylabel('|a|'); hold on;
    figure(102); subplot(2,1,1); semilogy(und_L_array,abs(Power_und(1:end-1))/1e3,'rs'); xlabel('Radiator s (m)'); ylabel('P_{out} (GW)'); hold on;
    figure(102); subplot(2,1,2); semilogy(und_L_array,abs(b_und(1:end-1)),'bs'); xlabel('Radiator s (m)'); ylabel('|b|'); hold on;
    figure(102); subplot(2,1,2); semilogy([0 und_L_array(end)],[1 1],'k--'); xlim([0 und_L_array(end)]); ylim([1e-5 10]); hold on;
end
%}
%{
% CALCULATE MBI-INDUCED SIDEBAND POWER
for m=1:1:length(Gf_11)
    %dnu_2rho(m)=(lambda_r-lambda_array(m)/C_ele(end)/100)/(2*rho_und*lambda_r);
    dnu_2rho(m)=lambda_r/(2*rho_und*lambda_array(m)/C_ele(end)/100);
    tmp15=1-sqrt(3)*z_sat+z_sat^2;
    tmp16=(2-2*sqrt(3)*z_sat+z_sat^2)/6;
    tmp17=(12-16*sqrt(3)*z_sat+22*z_sat^2-4*sqrt(3)*z_sat^3-z_sat^4)/108;
    if (dnu_2rho(m)<1)
        tmp18(m)=1;
    else
        tmp18(m)=0;
    end
end

Ps_d_fin_up=P_sat*A_hat_array.^2/9.*(tmp15-dnu_2rho*tmp16+dnu_2rho.^2*tmp17).*tmp18;
Ps_d_fin_lo=P_sat*A_hat_array.^2/9.*(tmp15+dnu_2rho*tmp16+dnu_2rho.^2*tmp17).*tmp18;
Ps_e_fin_up=P_sat*B_hat_array.^2/9*z_sat^2.*tmp18;
Ps_e_fin_lo=P_sat*B_hat_array.^2/9*z_sat^2.*tmp18;
Ps_fin_up=Ps_d_fin_up+Ps_e_fin_up;
Ps_fin_lo=Ps_d_fin_lo+Ps_e_fin_lo;

titletxt=sprintf('Main-signal power = %.2f MW',P_sat);

figure(104); stem(0,P_sat,'k'); hold on;
figure(104); stem(lambda_r./(lambda_array/C_ele(end)/100),Ps_fin_up,'b'); xlabel('\Delta\omega/\omega_0'); ylabel('Output power (MW)'); hold on;
figure(104); stem(-lambda_r./(lambda_array/C_ele(end)/100),Ps_fin_lo,'r'); title(titletxt); hold on;
%}

emit_gr=lambda_r/4/pi;

fprintf('========= SASE FEL PERFORMANCE ESTIMATE =========\n');
if (MX==0)
    fprintf('Based on 1-D model...\n');
else
    fprintf('Based on 3-D Ming-Xie model...\n');
end
fprintf('Saturation length z_sat = %.2f (%.2f) m \n',z_sat*Lg_und,und_lambda/rho_und);
fprintf('Saturation power P_sat = %.2f (%.2f) GW \n',P_sat,rho_und*egamma_vec(end)*0.511/1000*I_b_fin);
fprintf('Spectral width = %.2e \n',rho_und);
if (MX==0)
    fprintf('Power gain length = %.2f m \n',Lg_und);
else
    fprintf('Power gain length = %.2f (1-D %.2f) m \n',Lg_und,Lg_und_1D);
end
fprintf('Transverse radiation emittance = %.2f m \n',emit_gr);
fprintf('Transverse radiation size = %.2f m \n',sqrt(emit_gr*Lg_und));
fprintf('=================================================\n');

