function func=SES_calc(lambda_array,ekd_mat_fin_abs)

global I_b iIBS sigma_delta_IBS sigma_delta egamma C_ele

c_speed=2.99792458e10;        % unit: cm/sec
charge=1.6e-19;               % unit: C
mc2=0.51099895;               % unit: MeV

tmp01=2*charge*c_speed/I_b;
tmp02=egamma^3*mc2^2;
tmp03=trapz(lambda_array,(ekd_mat_fin_abs./lambda_array).^2);
tmp04=tmp01*tmp02*tmp03;      % unit: (MeV)^2,  SES due to collective kick
tmp05=(C_ele(end)*sigma_delta*egamma*mc2)^2;   %SES due to initial energy spread

p=0;
if (iIBS==1)
    while(isnan(sigma_delta_IBS(end-p)))
        p=p+1;
    end
    tmp06=(sigma_delta_IBS(end-p)*egamma*mc2)^2;%SES due to IBS
    func=sqrt(tmp05+tmp06+tmp04);
    %fprintf('SES contributed by IBS (%.1f%%) and by collective effect (%.1f%%)...\n',sqrt(tmp06)/func*100,sqrt(tmp04)/func*100);
    fprintf('%.4e\t%.4e\t%.4e keV...\n',sqrt(tmp05)*1e3,func*1e3,(sqrt(tmp06)-sqrt(tmp05))*1e3);
else
    func=sqrt(tmp05+tmp04);
    %fprintf('SES contributed by collective effect (%.1f%%)...\n',sqrt(tmp04)/func*100);
    fprintf('%.4e\t%.4e keV...\n',sqrt(tmp05)*1e3,func*1e3);
end
