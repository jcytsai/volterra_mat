function func=kernel_mod_w_TWLA(a,b)    
% this version modifies kernel function such that it can adopt the first
% argument a as a vector while keeping the second argument b as a scalar. 
% output of this subroutine: a vector func
%
% a: present time, a vector
% b: previous time, a scalar
%
    global k_wave s_ele C_ele sig_x_ele sig_y_ele betax0 betay0 sigma_delta re nb full_pipe_height;
    global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac egamma_vec emit_norm_x emit_norm_y;
    %global emitx emity egamma;
    global LSC_model ssCSR_model issCSRpp RF_info;
    
    format long
    
    egamma=interp1(s_ele,egamma_vec,a);  % vector
    emitx=emit_norm_x./egamma;
    emity=emit_norm_y./egamma;
    
    % base (x,px,y,py,z,delta)
    %tmp01=1i*k_wave*re*nb./egamma;    % constant term
    
    % base (x,px,y,py,z,Delta_gamma)
    tmp01=1i*k_wave*re*nb;
    
    tmp02=interp1(s_ele,C_ele,a);    % C(s), vector
    tmp03=interp1(s_ele,C_ele,b);    % C(s')
    
    % CSR-related models
    
    %{
    % Note: csr1d(k_wave*tmp03,b) <=> csr1d_nur(k_wave*tmp03,b,egamma)
    if ((iCSR_ss == 0)&&(iCSR_tr == 0)&&(iCSR_drift == 0))
        tmp04=0.0;
    elseif ((iCSR_ss == 1)&&(iCSR_tr == 0)&&(iCSR_drift == 0))
        tmp04=csr1d_nur(k_wave*tmp03,b,egamma);
    elseif ((iCSR_ss == 0)&&(iCSR_tr == 1)&&(iCSR_drift == 0))
        tmp04=csr1d_tr(k_wave*tmp03,b);
    elseif ((iCSR_ss == 0)&&(iCSR_tr == 0)&&(iCSR_drift == 1))
        tmp04=csr1d_drift(k_wave*tmp03,b);
    elseif ((iCSR_ss == 1)&&(iCSR_tr == 1)&&(iCSR_drift == 0))
        tmp04=csr1d_nur(k_wave*tmp03,b,egamma)+csr1d_tr(k_wave*tmp03,b);
    elseif ((iCSR_ss == 1)&&(iCSR_tr == 0)&&(iCSR_drift == 1))
        tmp04=csr1d_nur(k_wave*tmp03,b,egamma)+csr1d_drift(k_wave*tmp03,b);
    elseif ((iCSR_ss == 0)&&(iCSR_tr == 1)&&(iCSR_drift == 1))
        tmp04=csr1d_tr(k_wave*tmp03,b)+csr1d_drift(k_wave*tmp03,b);
    elseif ((iCSR_ss == 1)&&(iCSR_tr == 1)&&(iCSR_drift == 1))
        tmp04=csr1d_nur(k_wave*tmp03,b,egamma)+csr1d_tr(k_wave*tmp03,b)+csr1d_drift(k_wave*tmp03,b);
    end
    %}
    
    if (ssCSR_model==1); iNUR=0; iUR=1; end
    if (ssCSR_model==2); iNUR=1; iUR=0; end
    if (ssCSR_model==-1);iNUR=0; iUR=0; end
    
    if (iNUR==1 && iUR==1); fprintf('warning: both NUR- and UR-CSR models are used...\n'); end
    
    if (iCSR_ss==1)
        if (iNUR==1); 
            %tmp_CSR_ss_nur=csr1d_nur(k_wave*tmp03,b,egamma_b,iCSR_tr);
            tmp_CSR_ss_nur=csr1d_nur(k_wave*tmp03,b,egamma,iCSR_tr);
        end
        if (iUR==1);
            if (issCSRpp==1)
            	fac=1; 
                h=full_pipe_height;                    % full pipe height in cm (parameter added in GUI)
                tmp_rho=abs(auxr(b));
                if (tmp_rho<1e8)                       % within dipole
                    k_th=2*pi/fac*sqrt(tmp_rho/h^3);
                    if (k_wave*tmp03 < k_th)
                        if (first==1); fprintf('shielding CSR impedance calculated...\n'); first=first+1; end
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
    
    if (iCSR_tr==1)
        tmp_CSR_tr=csr1d_tr(k_wave*tmp03,b);
    else
        tmp_CSR_tr=0.0;
    end
    
    if (iCSR_drift==1)
        tmp_CSR_drift=csr1d_drift(k_wave*tmp03,b);
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
        a=3.5;                      % average iris radius in cm
        L=11.54;                    % cell/period length in cm
        g=9.0;                      % gap distance between irises in cm
        alpha=0.46;                 % scaling factor
        %}
        
        % CEBAF 7-cell SRF cavities (JLAB-TN-01-015)
        a=3.07;
        L=10.0;
        g=8.0;
        alpha=0.528;
        
        %{
        % SLAC 3-m TWLA
        a=11.6;
        L=3.5;
        g=2.9;
        alpha=0.5;
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
        if (mod(tmp11,2)==0)
            if (iLSC==1)
                tmp_LSC=lsc1d(k_wave*tmp03,interp1(s_ele,sig_x_ele,b),interp1(s_ele,sig_y_ele,b),b,LSC_model);
            else
                tmp_LSC=0.0;
            end
            tmp_LINAC=0.0;
        else
            % CEBAF 7-cell SRF cavities (JLAB-TN-01-015)
            a=3.07;
            L=10.0;
            g=8.0;
            alpha=0.528;
            tmp_LINAC=LINAC_imp(k_wave*tmp03,a,L,g,alpha,b);
            tmp_LSC=0.0;
        end
    end
    
    tmp04=tmp04+tmp_LSC+tmp_LINAC;
    %tmp04=csr1d_nur(k_wave*tmp03,b,egamma)+lsc1d(k_wave*tmp03,interp1(s_ele,sig_x_ele,b),interp1(s_ele,sig_y_ele,b));
    
    tmp05=R56transport(b,a);                                                      % R56(s'->s), vector
    
    tmp06=0.5*k_wave^2*emitx/(betax0);
    
    tmp07=betax0^2*R51mod(a,b).^2+R52mod(a,b).^2;                                 % vector
    %tmp07=betax0^2*R51mod_Deltagamma(a,b).^2+R52mod_Deltagamma(a,b).^2;          % vector
    
    %tmp08=0.5*(k_wave^2)*(sigma_delta^2)*(R56mod(a,b).^2);                       % vector
    %tmp08=0.5*(k_wave^2)*((sigma_delta*egamma).^2).*(R56mod(a,b).^2);            % vector
    tmp08=0.5*(k_wave^2)*((sigma_delta).^2).*(R56mod_Deltagamma(a,b).^2);         % vector
    
    tmp09=0.5*k_wave^2*emity/(betay0);
    
    tmp10=betay0^2*R53mod(a,b).^2+R54mod(a,b).^2;                                 % vector
    
    func=tmp01.*tmp02.*tmp03.*tmp04.*tmp05.*exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08);
    
    
