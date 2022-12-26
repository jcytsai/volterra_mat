    use_Piwinski=0;                      % Piwinski model
    use_CIMP=1;                          % Bjorken-Mtingwa model
    mod_factor=2*sqrt(2*log(2));         % for Gaussian bunch
    %mod_factor=1;                       % for uniform bunch
    update_SES=1;                        % update slice energy spread?
    
    scale_factor=1;
    tmp01=interp1(s_ele,C_ele,s);        % C(s)
    ds=s(2)-s(1);
    sigma_delta_IBS=zeros(1,length(s));
    %int_IBS_growth_rate_01=zeros(1,length(s));
    
    %clear tau_inv_IBS;
    if (use_CIMP==1)
        const_factor_01=re^2/(c_speed*charge);
        
        for m=1:1:length(s)
            
            if (m==1)
                sigma_delta_IBS(1,1)=sigma_delta;
            else
                %sigma_delta_IBS(1,m)=sigma_delta_IBS(1,1)*exp(int_IBS_growth_rate_01(1,m));
                sigma_delta_IBS(1,m)=sigma_delta_IBS(1,m-1)*exp(tau_inv_IBS(m-1)*ds);
            end
            
            s_loc=s(m);
            egamma_tmp=interp1(s_ele,egamma_vec,s_loc);
            betax=interp1(s_ele_Twiss,betax_s_ele,s_loc);
            betay=interp1(s_ele_Twiss,betay_s_ele,s_loc);
            alphax=interp1(s_ele_Twiss,alphax_s_ele,s_loc);
            alphay=interp1(s_ele_Twiss,alphay_s_ele,s_loc);
            R16=interp1(s_ele,R16_ele,s_loc);
            R26=interp1(s_ele,R26_ele,s_loc);
            R36=interp1(s_ele,R36_ele,s_loc);
            R46=interp1(s_ele,R46_ele,s_loc);
            sigx=interp1(s_ele,sig_x_ele,s_loc,'spline','extrap');
            sigy=interp1(s_ele,sig_y_ele,s_loc,'spline','extrap');
            emit_nx=interp1(s_ele,enx_s_ele,s_loc,'spline','extrap'); emit_gx=emit_nx/egamma_tmp;
            emit_ny=interp1(s_ele,eny_s_ele,s_loc,'spline','extrap'); emit_gy=emit_ny/egamma_tmp;
            
            %tmp12=2*mod_factor/(64*pi^2);  % 2: for single-pass accelerator (due to no synchrotron motion)
            tmp12=2/mod_factor/(64*pi^2);
            
            if (update_SES==0)
                A_coeff(m)=tmp12*const_factor_01*(I_b*tmp01(m))/(egamma_tmp^2*emit_nx*emit_ny*sigma_delta);
                %A_coeff(m)=tmp12*re^2*(nb*tmp01(m))/(egamma^2*emit_nx*emit_ny*sigma_delta); % same as above
            else
                A_coeff(m)=tmp12*const_factor_01*(I_b*tmp01(m))/(egamma_tmp^2*emit_nx*emit_ny*sigma_delta_IBS(1,m));
            end
            
            log_factor(m)=log(egamma_tmp^2*sigy*emit_gx/(re*betax));
            H_x(m)=(R16^2+(betax*R26+alphax*R16)^2)/betax;
            H_y(m)=(R36^2+(betay*R46+alphay*R36)^2)/betay;
            
            if (update_SES==0)
                sigma_H(m)=sqrt(1/(1/sigma_delta^2+H_x(m)/emit_gx+H_y(m)/emit_gy));
            else
                sigma_H(m)=sqrt(1/(1/sigma_delta_IBS(1,m)^2+H_x(m)/emit_gx+H_y(m)/emit_gy));
            end
            
            a_param(m)=sigma_H(m)/egamma_tmp*sqrt(betax/emit_gx);
            b_param(m)=sigma_H(m)/egamma_tmp*sqrt(betay/emit_gy);
            
            if (update_SES==0)
                tau_inv_IBS(m)=2*pi^(3/2)*A_coeff(m)*log_factor(m)*(sigma_H(m)^2/sigma_delta^2)*(Piwinski_g_func(b_param(m)/a_param(m))/a_param(m)+Piwinski_g_func(a_param(m)/b_param(m))/b_param(m));
            else
                tau_inv_IBS(m)=2*pi^(3/2)*A_coeff(m)*log_factor(m)*(sigma_H(m)^2/sigma_delta_IBS(1,m)^2)*(Piwinski_g_func(b_param(m)/a_param(m))/a_param(m)+Piwinski_g_func(a_param(m)/b_param(m))/b_param(m));
            end
            %{
            if (m<length(s) && m>=2)
                X=linspace(0,s(m),m);
                Y_01=tau_inv_IBS(1:m);
                int_IBS_growth_rate_01(1,m+1)=trapz(X,Y_01);
            end
            %}
        end
        %{
        if (update_SES==0)
            tau_inv_IBS=2*pi^(3/2)*A_coeff.*log_factor.*(sigma_H.^2/sigma_delta^2).*(Piwinski_g_func(b_param./a_param)./a_param+Piwinski_g_func(a_param./b_param)./b_param);
        else
            tau_inv_IBS=2*pi^(3/2)*A_coeff.*log_factor.*(sigma_H.^2./sigma_delta_IBS.^2).*(Piwinski_g_func(b_param./a_param)./a_param+Piwinski_g_func(a_param./b_param)./b_param);
        end
        %}
        %figure(999); plot(s/100,tau_inv_IBS*100); xlabel('s (m)'); ylabel('IBS growth rate (m^{-1})'); hold on;
        %figure(999); plot(s/100,H_x,'r'); xlabel('s (m)'); ylabel('H_x(red),H_y(blue) (cm)'); hold on;
        %figure(999); plot(s/100,H_y,'b'); xlabel('s (m)'); hold on;
        %figure(999); plot(s/100,log_factor,'r'); xlabel('s (m)'); ylabel('log factor'); hold on;
    end
    
    if (use_Piwinski==1)
        const_factor_01=re^2/(c_speed*charge);
        theta_max=0.01;         % 10 mrad, argument by SDM
        for m=1:1:length(s)
            
            if (m==1)
                sigma_delta_IBS(1,1)=sigma_delta;
            else
                %sigma_delta_IBS(1,m)=sigma_delta_IBS(1,1)*exp(int_IBS_growth_rate_01(1,m));
                sigma_delta_IBS(1,m)=sigma_delta_IBS(1,m-1)*exp(tau_inv_IBS(m-1)*ds);
            end
            
            s_loc=s(m);
            egamma_tmp=interp1(s_ele,egamma_vec,s_loc);
            R16=interp1(s_ele,R16_ele,s_loc);
            R26=interp1(s_ele,R26_ele,s_loc);
            R36=interp1(s_ele,R36_ele,s_loc);
            R46=interp1(s_ele,R46_ele,s_loc);
            
            if (R16==0)
                sigma_H(m)=sigma_delta_IBS(1,m);
                sigxp(m)=interp1(s_ele,sig_xp_ele,s_loc);
                if (sigxp(m)==0); sigxp(m)=min(abs(sig_xp_ele)); end
                emit_nx(m)=interp1(s_ele,enx_s_ele,s_loc,'spline','extrap');
                emit_ny(m)=interp1(s_ele,eny_s_ele,s_loc,'spline','extrap');
                if (emit_nx(1)==0); emit_nx(1)=enx_s_ele(1); end
                if (emit_ny(1)==0); emit_ny(1)=eny_s_ele(1); end
            else
                betax=interp1(s_ele_Twiss,betax_s_ele,s_loc);
                betay=interp1(s_ele_Twiss,betay_s_ele,s_loc);
                alphax=interp1(s_ele_Twiss,alphax_s_ele,s_loc);
                alphay=interp1(s_ele_Twiss,alphay_s_ele,s_loc);
                sigxp(m)=interp1(s_ele,sig_xp_ele,s_loc);
                %sigyp(m)=interp1(s_ele,sig_yp_ele,s_loc);
                emit_nx(m)=interp1(s_ele,enx_s_ele,s_loc,'spline','extrap');
                emit_ny(m)=interp1(s_ele,eny_s_ele,s_loc,'spline','extrap');
                
                if (sigxp(m)==0); sigxp(m)=min(abs(sig_xp_ele)); end
                if (emit_nx(1)==0); emit_nx(1)=enx_s_ele(1); end
                if (emit_ny(1)==0); emit_ny(1)=eny_s_ele(1); end
                
                emit_gx(m)=emit_nx(m)/egamma_tmp; 
                emit_gy(m)=emit_ny(m)/egamma_tmp;
                H_x(m)=(R16^2+(betax*R26+alphax*R16)^2)/betax;
                H_y(m)=(R36^2+(betay*R46+alphay*R36)^2)/betay;
                sigma_H(m)=sqrt(1/(1/sigma_delta_IBS(1,m)^2+H_x(m)/emit_gx(m)+H_y(m)/emit_gy(m)));
            end
            %tmp13=mod_factor/8;
            tmp13=1/(mod_factor*8);
            A_coeff(m)=tmp13*const_factor_01*(I_b*tmp01(m))/(egamma_tmp^2*emit_nx(m)*emit_ny(m)*sigma_delta_IBS(1,m)^2);
            log_factor(m)=log(egamma_tmp*sigxp(m)*emit_nx(m)*theta_max/(4*re));
            tau_inv_IBS(m)=A_coeff(m)*log_factor(m)*sigma_H(m)/sigma_delta_IBS(1,m);
        end
        %figure(999); plot(s/100,tau_inv_IBS*100); xlabel('s (m)'); ylabel('IBS growth rate (m^{-1})'); hold on;
    end
    %{
    prog=0; fprintf(1,'calculate integrated IBS growth: %1d%%\n',prog);
    for m=1:1:length(s)
    	if (m==1)
        	int_IBS_growth_rate_01(1)=0;
            %int_IBS_growth_rate_02(1)=0;
            %sigma_delta_IBS(1)=sigma_delta;
        else
            X=linspace(0,s(m),m);
            %Y_01=tau_inv_IBS(1:m);
            Y_01=tau_inv_IBS(1:m);
            int_IBS_growth_rate_01(m)=trapz(X,Y_01);
            %int_IBS_growth_rate_01(m)=tau_inv_IBS(m)*ds;
            %sigma_delta_IBS(m)=sigma_delta_IBS(1)*exp(int_IBS_growth_rate_01(m));
            %{
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
            else
                int_IBS_growth_rate_02=zeros(1,length(s));
            end
            %}
        end
        prog=100*m/length(s);
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    fprintf('\n');
    
    %int_IBS_growth_rate_01=scale_factor*int_IBS_growth_rate_01;
    %int_IBS_growth_rate_02=scale_factor*int_IBS_growth_rate_02;
    %}
    %figure(999); plot(s/100,int_IBS_growth_rate_01,'r'); xlabel('s (m)'); ylabel('integrated IBS growth (unitless)'); hold on;
    %figure(999); plot(s/100,int_IBS_growth_rate_02,'g'); xlabel('s (m)'); hold on;
    %figure(999); plot(s/100,sigma_delta_IBS,'g'); xlabel('s (m)'); ylabel('\sigma_{\delta} with IBS'); hold on;
    
    