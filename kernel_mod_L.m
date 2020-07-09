function func=kernel_mod_L(a,b)    

    global egamma k_wave s_ele C_ele sig_x_ele sig_y_ele emitx emity alphax0 alphay0 betax0 betay0 sigma_delta re nb full_pipe_height;
    global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac egamma_vec emit_norm_x emit_norm_y;
    global LSC_model ssCSR_model issCSRpp first RF_info;
    
    format long
    
    %ctrl_ss=dipole_action_control(b,2,iCSR_ss);
    %ctrl_tr=dipole_action_control(b,2,iCSR_tr);
    
    tmp01=re*nb./egamma;          % constant term
    %tmp02=interp1(s_ele,C_ele,a);             % C(s), vector
    tmp03=interp1(s_ele,C_ele,b);             % C(s')
    
    % CSR-related models
    
    if (ssCSR_model==1); iNUR=0; iUR=1; end   % ultra-relativistic CSR
    if (ssCSR_model==2); iNUR=1; iUR=0; end   % non-ultra-relativistic CSR
    if (ssCSR_model==-1);iNUR=0; iUR=0; end
    
    if (iNUR==1 && iUR==1) 
        fprintf('Warning: either NUR- or UR-CSR models can be used, not both...\n'); 
        pause; 
    end
    
    %if (iCSR_ss==1 && ctrl_ss==1)
    if (iCSR_ss==1)
        if (iNUR==1)
            tmp_CSR_ss_nur=csr1d_nur(k_wave*tmp03,b,egamma,iCSR_tr);
        end
        if (iUR==1)
            if (issCSRpp==1)
            	fac=1; 
                h=full_pipe_height;           % full pipe height in cm (parameter added in GUI)
                tmp_rho=abs(auxr(b));
                if (tmp_rho<1e8)              % within dipole
                    k_th=2*pi/fac*sqrt(tmp_rho/h^3);
                    if (k_wave*tmp03 < k_th)
                        if (first==1) 
                            fprintf('shielding CSR impedance used...\n'); 
                            first=first+1; 
                        end
                        tmp_CSR_ss_ur=csr_shielding_pp(k_wave,h,b);
                    else
                        tmp_CSR_ss_ur=csr1d(k_wave*tmp03,b,iCSR_tr);
                    end
                end
            else
                tmp_CSR_ss_ur=csr1d(k_wave*tmp03,b,iCSR_tr);
            end
        end
    else
        if (iNUR==1); tmp_CSR_ss_nur=0.0; end
        if (iUR==1);  tmp_CSR_ss_ur=0.0;  end
    end
    
    %if (iCSR_tr==1 && ctrl_tr==1)
    if (iCSR_tr==1)
        tmp_CSR_tr=csr1d_tr(k_wave*tmp03,b);
    else
        tmp_CSR_tr=0.0;
    end
    
    if (iCSR_drift==1)
        tmp_CSR_drift=csr1d_drift(k_wave*tmp03,b);
        %tmp_CSR_drift=csr1d_nur_drift(k_wave*tmp03,b);
    else
        tmp_CSR_drift=0.0;
    end
    
    if (iNUR==1); tmp04=tmp_CSR_ss_nur+tmp_CSR_tr+tmp_CSR_drift; end
    if (iUR==1);  tmp04=tmp_CSR_ss_ur+tmp_CSR_tr+tmp_CSR_drift; end
    if (iNUR==0 && iUR==0); tmp04=0.0; end
    
    tmp_LSC=0.0;
    tmp_LINAC=0.0;
    
    %{ 
    % when un-comment, remember to modify LINAC_imp.m
    % LSC-related models
    if (iLSC==1)
        tmp_LSC=lsc1d(k_wave*tmp03,interp1(s_ele,sig_x_ele,b),interp1(s_ele,sig_y_ele,b),b,LSC_model);
    else
        tmp_LSC=0.0;
    end
    
    
    % LINAC geometric impedance
    if (ilinac==1)
        %{
        % TESLA 9-cell SRF cavities
        a=3.5;                                % average iris radius in cm
        L=11.54;                              % cell/period length in cm
        g=9.0;                                % gap distance between irises in cm
        %alpha=0.46;                          % scaling factor
	alpha=1-0.465*sqrt(g/L)-0.07*(g/L);
        %}
        
        % CEBAF 7-cell SRF cavities (JLAB-TN-01-015)
        a=3.07;
        L=10.0;
        g=8.0;
        %alpha=0.528;
	alpha=1-0.465*sqrt(g/L)-0.07*(g/L);        

        %{
        % SLAC 3-m TWLA
        a=11.6;
        L=3.5;
        g=2.9;
        %alpha=0.5;
	alpha=1-0.465*sqrt(g/L)-0.07*(g/L);
        %}
        tmp_LINAC=LINAC_imp(k_wave*tmp03,a,L,g,alpha,b);
    else
        tmp_LINAC=0.0;
    end
    %}
    
    if ((iLSC==1)||(ilinac==1))
        if (ilinac==1)
            tmp11=find((RF_info(:,1)-b/100)<0,1,'last');
        else
            tmp11=0.0;
        end
        if (mod(tmp11,2)==0)                  % outside RFCA/TWLA
            if (iLSC==1)
                if (LSC_model==-1); fprintf('error: LSC model not specified...\n'); end
                tmp_LSC=lsc1d(k_wave*tmp03,interp1(s_ele,sig_x_ele,b),interp1(s_ele,sig_y_ele,b),b,LSC_model);
            else
                tmp_LSC=0.0;
            end
            tmp_LINAC=0.0;
        else                                  % inside RFCA/TWLA
            % CEBAF 7-cell SRF cavities (JLAB-TN-01-015)
            a=3.07;
            L=10.0;
            g=8.0;
            alpha=0.528;
            tmp_LINAC=LINAC_imp(k_wave*tmp03,a,L,g,alpha,b);
            %tmp_LINAC=0.0;
            %tmp_LSC=lsc1d(k_wave*tmp03,interp1(s_ele,sig_x_ele,b),interp1(s_ele,sig_y_ele,b),b,LSC_model);
            tmp_LSC=0.0;
        end
    end
    
    tmp04=tmp04+tmp_LSC+tmp_LINAC;                                             % Z(k0*C(s'),s'),scalar
    %tmp05=R56mod(a,b);                                                        % U(s,s')=R56(s)-R56(s'), vector
    tmp06=0.5*k_wave^2*emitx/(betax0);
    %tmp07=betax0^2*R51mod(a,b).^2+R52mod(a,b).^2;                             % vector
    tmp07=betax0^2*(R51mod(a,b)-R52mod(a,b)*alphax0/betax0).^2+R52mod(a,b).^2; % vector
    tmp08=0.5*(k_wave^2)*(sigma_delta^2)*(R56mod(a,b).^2);                     % vector
    tmp09=0.5*k_wave^2*emity/(betay0);
    %tmp10=betay0^2*R53mod(a,b).^2+R54mod(a,b).^2;                             % vector
    tmp10=betay0^2*(R53mod(a,b)-R54mod(a,b)*alphay0/betay0).^2+R54mod(a,b).^2; % vector
    
    func=tmp01.*tmp03.*tmp04.*exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08);
    
    
