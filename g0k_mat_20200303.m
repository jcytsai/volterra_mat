function func=g0k_mat(s)
% s: array from s_start to s_end
    global s_ele sig_x_ele sig_y_ele sig_xp_ele sig_yp_ele enx_s_ele eny_s_ele R16_ele R26_ele R36_ele R46_ele R51_ele R52_ele R53_ele R54_ele R56_ele;
    global C_ele I_b nb;
    global k_wave;
    global alphax0;
    global alphay0;
    global betax0;
    global betay0;
    global sigma_delta;
    global n_1k0;
    global emitx emity;
    %global egamma_vec egamma emit_norm_x emit_norm_y find_RF_para;
    
    global iIBS egamma s_ele_Twiss alphax_s_ele alphay_s_ele betax_s_ele betay_s_ele;
    
    format long
    
    re=2.81794032272e-13;                       % classical electron radius, cm
    c_speed=2.99792458e10;                      % speed of light, cm/sec
    charge=1.6e-19;
    
    ds=s(2)-s(1);
    
    %{
    %global emitx;
    %global emity;
    
    tmp_egamma_vec=interp1(s_ele,egamma_vec,s);
    emitx=emit_norm_x./tmp_egamma_vec;
    emity=emit_norm_y./tmp_egamma_vec;
    
    tmp01=interp1(s_ele,C_ele,s);               % C(s)
    tmp02=interp1(s_ele,R51_ele,s);             % R51(s)
    tmp03=interp1(s_ele,R52_ele,s);             % R52(s)
    tmp04=interp1(s_ele,R53_ele,s);             % R53(s)
    tmp05=interp1(s_ele,R54_ele,s);             % R54(s)
    tmp06=interp1(s_ele,R56_ele,s);             % R56(s)
    
    tmp07=(tmp01.^2)*(k_wave^2).*(emitx)/(2*betax0);
    tmp08=(betax0^2)*(tmp02.^2)+(tmp03.^2);
    tmp09=(tmp01.^2)*(k_wave^2).*(emity)/(2*betay0);
    tmp10=(betay0^2)*(tmp04.^2)+(tmp05.^2);
    
    
    if (find_RF_para==1)
        tmp11=(tmp01.^2).*(k_wave^2).*((sigma_delta*tmp_egamma_vec).^2).*(tmp06.^2)/2;
        %tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    else
        tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    end
    
    func=n_1k0*exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    %}
    
    %---------------------------------------------------------------------%
    % use initial egamma, PRST-AB 10, 104401 (2007)
    
    %tmp_egamma_vec=interp1(s_ele,egamma_vec,s);
    
    tmp01=interp1(s_ele,C_ele,s);               % C(s)
    tmp02=interp1(s_ele,R51_ele,s);             % R51(s)
    tmp03=interp1(s_ele,R52_ele,s);             % R52(s)
    tmp04=interp1(s_ele,R53_ele,s);             % R53(s)
    tmp05=interp1(s_ele,R54_ele,s);             % R54(s)
    tmp06=interp1(s_ele,R56_ele,s);             % R56(s)
    
    tmp07=(tmp01.^2)*(k_wave^2).*(emitx)/(2*betax0);
    %tmp08=(betax0^2)*(tmp02.^2)+(tmp03.^2);
    tmp08=(betax0^2)*((tmp02-tmp03*alphax0/betax0).^2)+(tmp03.^2);
    tmp09=(tmp01.^2)*(k_wave^2).*(emity)/(2*betay0);
    %tmp10=(betay0^2)*(tmp04.^2)+(tmp05.^2);
    tmp10=(betay0^2)*((tmp04-tmp05*alphay0/betay0).^2)+(tmp05.^2);
    tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    
    LD=exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    
    func=n_1k0*LD;
    %---------------------------------------------------------------------%
    
    use_Piwinski=0;                      % Piwinski model
    use_CIMP=0;                          % Bjorken-Mtingwa model
    HOM_IBS=0;                           % if 1, include second integral term, usually small
    mod_factor=2*sqrt(2*log(2));         % for Gaussian bunch
    %mod_factor=1;                        % for uniform bunch
    
    if (use_Piwinski==1 && use_CIMP==1)
        fprintf('cannot simultaeneously set use_Piwinski and use_CIMP, check g0k_mat...\n');
    end
    
    if (iIBS==1 && use_CIMP==1)
        const_factor_01=re^2/(c_speed*charge);
        
        for m=1:1:length(s)
            s_loc=s(m);
            betax=interp1(s_ele_Twiss,betax_s_ele,s_loc);
            betay=interp1(s_ele_Twiss,betay_s_ele,s_loc);
            alphax=interp1(s_ele_Twiss,alphax_s_ele,s_loc);
            alphay=interp1(s_ele_Twiss,alphay_s_ele,s_loc);
            R16=interp1(s_ele,R16_ele,s_loc);
            R26=interp1(s_ele,R26_ele,s_loc);
            R36=interp1(s_ele,R36_ele,s_loc);
            R46=interp1(s_ele,R46_ele,s_loc);
            sigx=interp1(s_ele,sig_x_ele,s_loc);
            sigy=interp1(s_ele,sig_y_ele,s_loc);
            emit_nx=interp1(s_ele,enx_s_ele,s_loc); emit_gx=emit_nx/egamma;
            emit_ny=interp1(s_ele,eny_s_ele,s_loc); emit_gy=emit_ny/egamma;
            
            tmp12=2*mod_factor/(64*pi^2);  % 2: for single-pass accelerator (due to no synchrotron motion)
            A_coeff(m)=tmp12*const_factor_01*(I_b*tmp01(m))/(egamma^2*emit_nx*emit_ny*sigma_delta);
            %A_coeff(m)=tmp12*re^2*(nb*tmp01(m))/(egamma^2*emit_nx*emit_ny*sigma_delta); % same as above
            log_factor(m)=log(egamma^2*sigy*emit_gx/(re*betax));
            H_x=(R16^2+(betax*R26+alphax*R16)^2)/betax;
            H_y=(R36^2+(betay*R46+alphay*R36)^2)/betay;
            sigma_H(m)=sqrt(1/(1/sigma_delta^2+H_x/emit_gx+H_y/emit_gy));
            a_param(m)=sigma_H(m)/egamma*sqrt(betax/emit_gx);
            b_param(m)=sigma_H(m)/egamma*sqrt(betay/emit_gy);
        end
        
        tau_inv_IBS=2*pi^(3/2)*A_coeff.*log_factor.*(sigma_H.^2/sigma_delta^2).*(Piwinski_g_func(b_param./a_param)./a_param+Piwinski_g_func(a_param./b_param)./b_param);
        %figure(999); plot(s/100,tau_inv_IBS*100); xlabel('s (m)'); ylabel('IBS growth rate (m^{-1})'); hold on;
        
    end
    
    if (iIBS==1 && use_Piwinski==1)
        const_factor_01=re^2/(c_speed*charge);
        theta_max=0.01;         % 10 mrad, argument by SDM
        for m=1:1:length(s)
            s_loc=s(m);
            R16=interp1(s_ele,R16_ele,s_loc);
            R26=interp1(s_ele,R26_ele,s_loc);
            R36=interp1(s_ele,R36_ele,s_loc);
            R46=interp1(s_ele,R46_ele,s_loc);
            
            if (R16==0)
                sigma_H(m)=sigma_delta;
                sigxp(m)=interp1(s_ele,sig_xp_ele,s_loc);
                if (sigxp(m)==0); sigxp(m)=min(abs(sig_xp_ele)); end
            else
                betax=interp1(s_ele_Twiss,betax_s_ele,s_loc);
                betay=interp1(s_ele_Twiss,betay_s_ele,s_loc);
                alphax=interp1(s_ele_Twiss,alphax_s_ele,s_loc);
                alphay=interp1(s_ele_Twiss,alphay_s_ele,s_loc);
                sigxp(m)=interp1(s_ele,sig_xp_ele,s_loc);
                %sigyp(m)=interp1(s_ele,sig_yp_ele,s_loc);
                emit_nx(m)=interp1(s_ele,enx_s_ele,s_loc);
                emit_ny(m)=interp1(s_ele,eny_s_ele,s_loc);
                
                if (sigxp(m)==0); sigxp(m)=min(abs(sig_xp_ele)); end
                if (emit_nx(1)==0); emit_nx(1)=enx_s_ele(1); end
                if (emit_ny(1)==0); emit_ny(1)=eny_s_ele(1); end
                
                emit_gx(m)=emit_nx(m)/egamma; 
                emit_gy(m)=emit_ny(m)/egamma;
                H_x=(R16^2+(betax*R26+alphax*R16)^2)/betax;
                H_y=(R36^2+(betay*R46+alphay*R36)^2)/betay;
                sigma_H(m)=sqrt(1/(1/sigma_delta^2+H_x/emit_gx(m)+H_y/emit_gy(m)));
            end
        end
        tmp13=mod_factor/8;
        A_coeff=tmp13*const_factor_01*(I_b*tmp01)./(egamma^2*emit_nx.*emit_ny*sigma_delta^2);
        log_factor=log(egamma*sigxp.*emit_nx*theta_max/(4*re));
        tau_inv_IBS=A_coeff.*log_factor.*sigma_H/sigma_delta;
        %figure(999); plot(s/100,tau_inv_IBS*100); xlabel('s (m)'); ylabel('IBS growth rate (m^{-1})'); hold on;
    end
    
    if (iIBS==1)
    for m=1:1:length(s)
    	if (m==1)
        	int_IBS_growth_rate_01(1)=0;
            int_IBS_growth_rate_02(1)=0;
        else
            X=linspace(0,s(m),m);
            %Y_01=tau_inv_IBS(1:m);
            Y_01=tau_inv_IBS(1:m);
            %int_IBS_growth_rate_01(m)=trapz(X,Y_01);
            int_IBS_growth_rate_01(m)=tau_inv_IBS(m)*ds;
            if (HOM_IBS==1)
            	for n=1:1:m
                    if (n==m)
                    	R56_rel(n)=0;
                    else
                        R56_rel(n)=R56transport(s(n),s(m));
                    end
                end
                Y_02=tau_inv_IBS(1:m).*R56_rel(1:m);  % a: s, vector; b: s', scalar
                %int_IBS_growth_rate_02(m)=trapz(X,Y_02);
                int_IBS_growth_rate_02(m)=tau_inv_IBS(m)*R56_rel(m)*ds;
            end
        end
    end
    %figure(999); plot(s/100,int_IBS_growth_rate_01,'r'); xlabel('s (m)'); ylabel('integrated IBS growth (unitless)'); hold on;
    %figure(999); plot(s/100,int_IBS_growth_rate_02,'g'); xlabel('s (m)'); hold on;
    
    if (HOM_IBS==1)
        LD_mod=(tmp01.^2)*(k_wave^2).*sigma_delta^2.*tmp06;
        func=n_1k0*LD.*(1-int_IBS_growth_rate_01+LD_mod.*int_IBS_growth_rate_02);
    else
        func=n_1k0*LD.*(1-int_IBS_growth_rate_01);
    end
    end
    
    