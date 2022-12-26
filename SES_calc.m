function [func_01,func_02,func_03]=SES_calc(lambda_array,ekd_mat_fin_abs,Gf_11,Gc_total)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: load a case, then
% run SES_calc(lambda_array,ekd_mat_fin_abs,Gf_11,Gc_total);     % unit: MeV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global I_b I_A sigma_delta egamma egamma_vec s_ele C_ele sig_x_ele sig_y_ele nb R56_ele emit_norm_x emit_norm_y;
global iIBS s_IBS sigma_delta_IBS iSES_ana chirp;

c_speed=2.99792458e10;        % unit: cm/sec
%charge=1.6e-19;               % unit: C
charge=4.803204e-10;          % unit: esu
mc2=0.51099895;               % unit: MeV
re=2.8179403267e-13;          % unit: cm

stage_num=2;                  % multistage number

egamma_vec_ext=interp1(s_ele,egamma_vec,s_IBS);

lambda_array_tmp=linspace(lambda_array(1),2*lambda_array(end),1000);
ekd_mat_fin_abs_tmp=interp1(lambda_array,ekd_mat_fin_abs,lambda_array_tmp,'spline','extrap');

% avoid data zero-crossing
[val,ind]=min(abs(ekd_mat_fin_abs_tmp(length(lambda_array):end)));
if (val<=0.0)
    lambda_array_ext=lambda_array_tmp(1:length(lambda_array)+ind-1);
    ekd_mat_fin_abs_ext=ekd_mat_fin_abs_tmp(1:length(lambda_array)+ind-1);
else
    lambda_array_ext=lambda_array_tmp;
    ekd_mat_fin_abs_ext=ekd_mat_fin_abs_tmp;
end



%{
tmp01=2*charge*c_speed/I_b;
tmp02=egamma^3*mc2^2;
tmp03=trapz(lambda_array,(ekd_mat_fin_abs./lambda_array).^2);
tmp04=tmp01*tmp02*tmp03;      % unit: (MeV)^2,  SES due to collective kick
%}
%tmp02=(egamma*mc2)^2;
tmp02=(mc2)^2;                % 20200409, after benchmark with G.Perosa's note
%tmp03=trapz(lambda_array,(ekd_mat_fin_abs./lambda_array).^2);
tmp03=trapz(lambda_array_ext,(ekd_mat_fin_abs_ext./lambda_array_ext).^2);

%tmp04=2/nb*tmp02*tmp03*(egamma_vec(end))^2;              % unit: (MeV)^2,  SES due to collective kick
%tmp04=2*charge*c_speed/I_b*tmp03;                        % 2*charge*c_speed/I_b = 2/nb, 20200905

%
mod_factor_01=1;             % correction of Gaussian energy width? maybe unnecessary, 20200920
mod_factor_02=4*pi/376.73;   % convert impedance from CGS to MKS
Z_LSC=zeros(length(lambda_array),length(s_IBS));
k_wave=2*pi./lambda_array;
C_s_IBS=interp1(s_ele,C_ele,s_IBS);
tmp07=I_b*C_s_IBS./egamma_vec_ext/I_A;
for m=1:1:length(lambda_array)
    Z_LSC(m,:)=lsc1d_vec(k_wave(m)*C_s_IBS,interp1(s_ele,sig_x_ele,s_IBS),interp1(s_ele,sig_y_ele,s_IBS),s_IBS,1);
    tmp08(m)=trapz(s_IBS,tmp07.*Z_LSC(m,:).*Gc_total(:,m).');
    tmp081(m)=trapz(s_IBS,C_s_IBS.*Z_LSC(m,:).*Gc_total(:,m).')/(tmp07(end)^stage_num);       % remove I_b/egamma/I_A
    tmp10(m)=trapz(s_IBS,abs(Z_LSC(m,:)));
end
[~,ind]=max(Gf_11);
%tmp09=trapz(lambda_array,abs(tmp08).^2./(lambda_array.^2));
tmp09=trapz(lambda_array(1:ind),abs(tmp08(1:ind)).^2./(lambda_array(1:ind).^2));
tmp11_int=abs(tmp081).^2./(lambda_array.^2);
%tmp04=8/nb*tmp09*tmp02*(egamma_vec(end))^2/mod_factor_01*mod_factor_02^2;
tmp04=8/nb*tmp09*tmp02*(egamma_vec(end))^2/2000*C_ele(end);
tmp091=trapz(lambda_array,abs(tmp08).^2./(lambda_array.^2));
tmp041=8/nb*tmp091*tmp02*(egamma_vec(end))^2/mod_factor_01*mod_factor_02^2;
%}

%{
% From Eq.(19), S. Di Mitri, NJP (2020), something incorrect
for m=1:1:length(lambda_array)
    Z_LSC(m,:)=lsc1d_vec(k_wave(m)*C_s_IBS,interp1(s_ele,sig_x_ele,s_IBS),interp1(s_ele,sig_y_ele,s_IBS),s_IBS,1);
    tmp08(m)=trapz(s_IBS,Z_LSC(m,:));
end
tmp09=trapz(lambda_array,abs(Gf_11.*tmp08).^2./(lambda_array.^2));
tmp04=2/nb*tmp09*tmp02*(egamma_vec(end))^2;
%}

tmp05=(C_ele(end)*sigma_delta*egamma*mc2)^2;             % unit: (MeV)^2,  SES due to initial energy spread

%figure(123); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array_ext*10^4,(ekd_mat_fin_abs_ext./lambda_array_ext).^2,'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|p(k_z;s_f)/\lambda|^2'); hold on;
figure(123); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,tmp11_int,'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('integrand'); hold on;
%figure(124); yyaxis left; plot(lambda_array*10^4,Gf_11); hold on;
%figure(124); yyaxis right; plot(lambda_array*10^4,tmp10); hold on;

% FIND THRESHOLD CURRENT FOR WHICH SES WITH IBS = SES WITHOUT IBS
% USE GRAPHICAL SOLUTION

C_tot=C_ele(end);
C_BC1=2;               % valid only for identical two-BC chicane

if (C_tot==1)
    factor=1.0;        % without bunch compression
else
    factor=0.75;       % valid only for identical two-BC chicane
end

[sigmaNew,muNew,Anew]=mygaussfit(lambda_array,tmp11_int);
y=Anew*exp(-(lambda_array-muNew).^2/(2*sigmaNew^2));
figure(123); plot(lambda_array*1e4,y,'.r'); hold on;
Anew=Anew/(2*pi);

I_b_array=linspace(0,50,1000);                    % unit: Amp
tmp12=2/(2*sqrt(2*log(2)))*re^2*s_IBS(end)*factor;
tmp13=c_speed*charge*mean(sig_x_ele)*emit_norm_x;
%tmp13=c_speed*charge*mean(sig_x_ele)*mean(sig_y_ele)*emit_norm_x*emit_norm_y;
tmp14=Anew*8/nb/2000/(I_A^(2+2*stage_num)*egamma^(2*stage_num))*sqrt(pi/2)*sigmaNew;       % a factor of 2000 for the same reason as tmp04
B_coeff=2*pi^2*(R56_ele(end))^2*tmp12/tmp13/(egamma_vec(end))^2*I_b_array*C_tot^2;
D_coeff=1+6*B_coeff*sigmaNew^2/muNew^4;
tmp15=muNew/(sqrt(2)*sigmaNew);
tmp16=1+erf(tmp15);
tmp17=1+erf(sqrt(D_coeff)*tmp15);
tmp18=sqrt(D_coeff).*exp(B_coeff/(muNew^2));

LHS=tmp12/tmp13*I_b_array;
RHS=tmp14*C_tot^(2*stage_num)*I_b_array.^(2+2*stage_num).*(tmp16-tmp17./tmp18);
[~,ind]=min(abs(LHS(200:end)-RHS(200:end)));   % MAKE SURE TO MODIFY THE OFFSET THE NEXT LINE
func_03=I_b_array(200+ind);

%
figure(125); semilogy(I_b_array,LHS,'b'); xlabel('I_b (A)'); hold on;
figure(125); semilogy(I_b_array,RHS,'r'); xlabel('I_b (A)'); hold on;
%}
p=0;
if (iIBS==1)
    while(isnan(sigma_delta_IBS(end-p)))
        p=p+1;
    end
    tmp06=(sigma_delta_IBS(end-p)*egamma*mc2)^2;  % unit: (MeV)^2, SES due to IBS
    %func=sqrt(tmp06+tmp04);
    func_01=sqrt(tmp06+(C_ele(end))^2*tmp04);
    func_02=sqrt(tmp06+(C_ele(end))^2*tmp041);
    fprintf('final pure optics SES %.4e keV, total SES %.4e (%.4e) keV, extra IBS SES %.4e keV...\n',sqrt(tmp05)*1e3,func_01*1e3,func_02*1e3,(sigma_delta_IBS(end)-sigma_delta)*egamma*mc2*1e3);
else
    %func=sqrt(tmp05+tmp04);
    func_01=sqrt(tmp05+(C_ele(end))^2*tmp04);
    func_02=sqrt(tmp05+(C_ele(end))^2*tmp041);
    fprintf('final pure optics SES %.4e keV, total SES %.4e (%.4e) keV...\n',sqrt(tmp05)*1e3,func_01*1e3,func_02*1e3);
end
