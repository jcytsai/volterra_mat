function func=kernel_mod_IBS_zz(a,b,kernel_order)
% a: s, vector; b: s' or tau, scalar
% output of this subroutine: a vector func
    global egamma k_wave s_ele C_ele sig_x_ele sig_y_ele emitx emity alphax0 alphay0 betax0 betay0 sigma_delta re nb full_pipe_height;
    global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac egamma_vec emit_norm_x emit_norm_y;
    global LSC_model ssCSR_model issCSRpp first RF_info;
    global iIBS s_IBS tau_inv_IBS tau_inv_IBS_x tau_inv_IBS_y emit_gx_IBS emit_gy_IBS D_zz sigma_delta_IBS;
    
    format long
    
    [~,ind]=min(abs(b-s_IBS));
    
    tmp01=R56transport(b,a);                                                                         % R56(s'->s), vector
    tmp02=interp1(s_ele,C_ele,a);        % C(s), vector
    tmp03=interp1(s_ele,C_ele,b);        % C(s') or C(\tau), scalar
    %tmp06=0.5*k_wave^2*emitx/(betax0);
    tmp06=0.5*k_wave^2*emit_gx_IBS(ind)/(betax0);
    tmp07=betax0^2*(R51transport(b,a)-R52transport(b,a)*alphax0/betax0).^2+R52transport(b,a).^2;     % vector
    %tmp08=0.5*(k_wave^2)*(sigma_delta^2)*(tmp01.^2);                                                 % vector
    tmp08=0.5*(k_wave^2)*(sigma_delta_IBS(ind)^2)*(tmp01.^2);
    %tmp09=0.5*k_wave^2*emity/(betay0);
    tmp09=0.5*k_wave^2*emit_gy_IBS(ind)/(betay0);
    tmp10=betay0^2*(R53transport(b,a)-R54transport(b,a)*alphay0/betay0).^2+R54transport(b,a).^2;     % vector
    %tmp11=exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08);
    tmp11=exp(tmp02.^2.*(-tmp06.*tmp07-tmp09.*tmp10-tmp08));
    tmp12=interp1(s_IBS,D_zz,b);
    
    switch kernel_order
        case 1
            %func=tmp12*k_wave*tmp01.*tmp11;
            func=tmp12*k_wave*tmp02.*tmp01.*tmp11;
        case 2
            %func=tmp12*k_wave^2.*(tmp01.^2).*tmp11;
            func=tmp12*k_wave^2*(tmp02.^2).*(tmp01.^2).*tmp11;
        case 30
            %func=tmp12*k_wave^3*(tmp01.^3).*tmp11;
            func=tmp12*k_wave^3*(tmp02.^3).*(tmp01.^3).*tmp11;
        case 31
            %func=interp1(s_IBS,D_zz,b)*k_wave^3*(tmp01.^3).*exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08)*(sigma_delta^2);
            %func=tmp12*k_wave^3*(tmp01.^3).*tmp11*(sigma_delta_IBS(ind)^2);
            func=tmp12*k_wave^3*(tmp02.^3).*(tmp01.^3).*tmp11*(sigma_delta_IBS(ind)^2);
        case 4
            %func=interp1(s_IBS,D_zz,b)*k_wave^4*(tmp01.^4).*exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08)*(sigma_delta^2);
            %func=tmp12*k_wave^4*(tmp01.^4).*tmp11*(sigma_delta_IBS(ind)^2);
            func=tmp12*k_wave^4*(tmp02.^4).*(tmp01.^4).*tmp11*(sigma_delta_IBS(ind)^2);
    end
    