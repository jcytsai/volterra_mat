function func_01=SES_calc(lambda_array,ekd_mat_fin_abs,Gf_11,Gc_total)

global I_b I_A sigma_delta egamma egamma_vec s_ele C_ele sig_x_ele sig_y_ele nb R56_ele;
global iIBS s_IBS sigma_delta_IBS iSES_ana;

c_speed=2.99792458e10;        % unit: cm/sec
charge=1.6e-19;               % unit: C
mc2=0.51099895;               % unit: MeV

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
    tmp10(m)=trapz(s_IBS,abs(Z_LSC(m,:)));
end
[~,ind]=max(Gf_11);
%tmp09=trapz(lambda_array,abs(tmp08).^2./(lambda_array.^2));
tmp09=trapz(lambda_array(1:ind),abs(tmp08(1:ind)).^2./(lambda_array(1:ind).^2));
%tmp04=8/nb*tmp09*tmp02*(egamma_vec(end))^2/mod_factor_01*mod_factor_02^2;
tmp04=8/nb*tmp09*tmp02*(egamma_vec(end))^2/2000*C_ele(end);
%tmp04=8/nb*tmp09*tmp02*(egamma_vec(end))^2;
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
figure(123); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array*10^4,abs(tmp08).^2./(lambda_array.^2),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('integrand'); hold on;
figure(124); yyaxis left; plot(lambda_array*10^4,Gf_11); hold on;
figure(124); yyaxis right; plot(lambda_array*10^4,tmp10); hold on;

p=0;
if (iIBS==1)
    while(isnan(sigma_delta_IBS(end-p)))
        p=p+1;
    end
    %tmp06=(sigma_delta_IBS(end-p)*egamma_vec(end)*mc2)^2;  % unit: (MeV)^2, SES due to IBS
    tmp06=(sigma_delta_IBS(end-p)*egamma*mc2)^2;  % unit: (MeV)^2, SES due to IBS, 20201004
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

