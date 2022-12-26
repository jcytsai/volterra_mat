function func=kernel_mod_IBS_z(a,b,kernel_order)
% a: s, vector; b: s' or tau, scalar
% output of this subroutine: a vector func
    global egamma k_wave s_ele C_ele sig_x_ele sig_y_ele emitx emity alphax0 alphay0 betax0 betay0 sigma_delta re nb full_pipe_height;
    global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac egamma_vec emit_norm_x emit_norm_y;
    global LSC_model ssCSR_model issCSRpp first RF_info;
    global iIBS s_IBS tau_inv_IBS tau_inv_IBS_x tau_inv_IBS_y emit_gx_IBS emit_gy_IBS D_IBS_z sigma_delta_IBS;
    
    format long
    
    N_max=0;      % 0, no sum
    [~,ind]=min(abs(b-s_IBS));
    
    tmp01=R56transport(b,a);                                                                         % R56(s'->s), vector
    tmp02=interp1(s_ele,C_ele,a);        % C(s), vector
    tmp03=interp1(s_ele,C_ele,b);        % C(s') or C(\tau), scalar
    %tmp06=0.5*k_wave^2*emitx/(betax0);
    tmp06=0.5*k_wave^2*emit_gx_IBS(ind)/(betax0);
    tmp07=betax0^2*(R51transport(b,a)-R52transport(b,a)*alphax0/betax0).^2+R52transport(b,a).^2;     % vector
    %tmp08=0.5*(k_wave^2)*(sigma_delta^2)*(tmp01.^2);                                                % vector
    tmp08=0.5*(k_wave^2)*(sigma_delta_IBS(ind)^2)*(tmp01.^2);
    %tmp09=0.5*k_wave^2*emity/(betay0);
    tmp09=0.5*k_wave^2*emit_gy_IBS(ind)/(betay0);
    tmp10=betay0^2*(R53transport(b,a)-R54transport(b,a)*alphay0/betay0).^2+R54transport(b,a).^2;     % vector
    %tmp11=exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08);                                                     % {L.D.;s'->s} in 6D
    tmp11=exp(tmp02.^2.*(-tmp06.*tmp07-tmp09.*tmp10-tmp08));
    %tmp12=exp(-tmp06.*tmp07-tmp09.*tmp10);                                                           % {L.D.;s'->s} in 4D
    tmp12=exp(tmp02.^2.*(-tmp06.*tmp07-tmp09.*tmp10));
    %tmp13=erfi(0.5*k_wave*tmp01*sigma_delta_IBS(ind));
    %tmp13=erfi(0.5*k_wave*tmp02.*tmp01*sigma_delta_IBS(ind));
    tmp13=2*(0.5*k_wave*tmp02.*tmp01*sigma_delta_IBS(ind))/sqrt(pi);                                 % erfi with small argument
    
    if (N_max>0)
        tmp14=sum_hpg(N_max,-tmp08,1.5);    % "sum+3/2", vector
        tmp15=sum_hpg(N_max,-tmp08,2.5);    % "sum+5/2", vector
    elseif (N_max==0)
        tmp14=gamma(3/2)*cos(sqrt(abs(tmp08)*3));
        tmp15=gamma(5/2)*sinc(sqrt(abs(tmp08)*5)/pi);
    end
    
    tmp16=interp1(s_IBS,D_IBS_z,b);
    
    switch kernel_order
        case 0
            func=tmp16*tmp11.*tmp13;
        case 1
            %func=tmp16*k_wave*tmp01.*tmp11.*tmp13;
            func=tmp16*k_wave*tmp02.*tmp01.*tmp11.*tmp13;
        case 0.5
            func=tmp16*tmp12*sqrt(2/pi)/sigma_delta_IBS(ind).*tmp14;
        case 1.5
            %func=tmp16*k_wave*tmp01.*tmp12*sqrt(2/pi)*sigma_delta_IBS(ind).*tmp14;
            func=tmp16*k_wave*tmp02.*tmp01.*tmp12*sqrt(2/pi)*sigma_delta_IBS(ind).*tmp14;
        case 2.5
            %func=tmp16*(k_wave*tmp01).^2.*tmp12*2*sqrt(2)/pi*sigma_delta_IBS(ind).*tmp15;
            func=tmp16*(k_wave*tmp02.*tmp01).^2.*tmp12*2*sqrt(2)/pi*sigma_delta_IBS(ind).*tmp15;
    end
    