function func=SES_calc(lambda_array,ekd_mat_fin_abs,Gf_11)

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
tmp04=2/nb*tmp02*tmp03*(egamma_vec(end))^2;              % unit: (MeV)^2,  SES due to collective kick
tmp05=(C_ele(end)*sigma_delta*egamma*mc2)^2;             % unit: (MeV)^2,  SES due to initial energy spread

figure(1);  set(gca,'FontSize',40,'linewidth',5); plot(lambda_array_ext*10^4,(ekd_mat_fin_abs_ext./lambda_array_ext).^2,'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|p(k_z;s_f)/\lambda|^2'); hold on;

%iSES_ana=1;
if (iSES_ana==1)
    Gf_11_tmp=interp1(lambda_array,Gf_11,lambda_array_ext,'spline','extrap');
    
    % avoid data zero-crossing
    lambda_array_tmp=lambda_array_ext;
    [val,ind]=min(abs(Gf_11_tmp(length(lambda_array):end)));
    if (val<=0.0)
        lambda_array_ext=lambda_array_tmp(1:length(lambda_array)+ind-1);
        Gf_11_ext=Gf_11_tmp(1:length(lambda_array)+ind-1);
    else
        lambda_array_ext=lambda_array_tmp;
        Gf_11_ext=Gf_11_tmp;
    end
    
    Gf_11_ext=Gf_11_ext.*(Gf_11_ext>=0);
    figure(9);  set(gca,'FontSize',40,'linewidth',5); plot(lambda_array_ext*10^4,Gf_11_ext,'b-','linewidth',5); hold on;
    
    fprintf('SES postprocess done...%3d%%\n',0);
    for m=1:1:length(lambda_array_ext)
        for n=1:1:length(s_IBS)
            tmp07=interp1(s_ele,C_ele,s_IBS(n));
            Z_LSC(n)=lsc1d(2*pi/lambda_array_ext(m)*tmp07,interp1(s_ele,sig_x_ele,s_IBS(n)),interp1(s_ele,sig_y_ele,s_IBS(n)),s_IBS(n),1);
            tmp08(n)=I_b/I_A*tmp07*abs(Z_LSC(n));
        end
        Z_int(m)=trapz(s_IBS,tmp08);
        Z_int_egamma(m)=trapz(s_IBS,tmp08./egamma_vec_ext);
        
        prog=(100*(m/length(lambda_array_ext)));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');
    tmp11=2*pi./lambda_array_ext*C_ele(end)*abs(R56_ele(end));
    Gf_11_ext_ana=tmp11.*Z_int_egamma.*exp(-0.5*(tmp11*sigma_delta).^2);
    tmp09=(abs(Gf_11_ext.*Z_int)./lambda_array_ext).^2;
    tmp12=(abs(Gf_11_ext_ana.*Z_int)./lambda_array_ext).^2;
    tmp10=2/nb*mc2^2*trapz(lambda_array_ext,tmp09);        % unit: (MeV)^2,  SES due to collective kick, semi-analytical Gf_11
    tmp13=2/nb*mc2^2*trapz(lambda_array_ext,tmp12);        % unit: (MeV)^2,  SES due to collective kick, analytical Gf_11
end

p=0;
if (iIBS==1)
    while(isnan(sigma_delta_IBS(end-p)))
        p=p+1;
    end
    tmp06=(sigma_delta_IBS(end-p)*egamma_vec(end)*mc2)^2;  % unit: (MeV)^2, SES due to IBS
    %func=sqrt(tmp05+tmp06+tmp04);
    func=sqrt(tmp06+tmp04);
    if (iSES_ana==1)
        func_ana=sqrt(tmp06+tmp10);
        func_ana_ana=sqrt(tmp06+tmp13);
    end
    if (iSES_ana==1)
        fprintf('pure optics SES %.4e keV, total SES %.4e keV (analytical %.4e keV, %.4e keV), extra IBS SES %.4e keV...\n',sqrt(tmp05)*1e3,func*1e3,func_ana*1e3,func_ana_ana*1e3,sigma_delta*egamma*mc2*1e3);
    else
        %fprintf('SES contributed by IBS (%.1f%%) and by collective effect (%.1f%%)...\n',sqrt(tmp06)/func*100,sqrt(tmp04)/func*100);
        fprintf('pure optics SES %.4e keV, total SES %.4e keV, extra IBS SES %.4e keV...\n',sqrt(tmp05)*1e3,func*1e3,sigma_delta*egamma*mc2*1e3);
    end
else
    func=sqrt(tmp05+tmp04);
    if (iSES_ana==1)
        func_ana=sqrt(tmp05+tmp10);
        func_ana_ana=sqrt(tmp05+tmp13);
    end
    if (iSES_ana==1)
        fprintf('pure optics SES %.4e keV, total SES %.4e keV (analytical %.4e keV, %.4e keV)...\n',sqrt(tmp05)*1e3,func*1e3,func_ana*1e3,func_ana_ana*1e3);
    else
        %fprintf('SES contributed by collective effect (%.1f%%)...\n',sqrt(tmp04)/func*100);
        fprintf('pure optics SES %.4e keV, total SES %.4e keV...\n',sqrt(tmp05)*1e3,func*1e3);
    end
end

