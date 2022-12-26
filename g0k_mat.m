function func=g0k_mat(s)

    global s_ele;
    global R51_ele;
    global R52_ele;
    global R53_ele;
    global R54_ele;
    global R56_ele;
    global C_ele;
    global k_wave;
    global alphax0;
    global alphay0;
    global betax0;
    global betay0;
    global sigma_delta;
    global n_1k0;
    %global egamma_vec egamma emit_norm_x emit_norm_y find_RF_para;
    global iIBS s_IBS emit_gx_IBS emit_gy_IBS sigma_delta_IBS;
    
    format long
    
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
    
    global emitx emity;
    
    tmp01=interp1(s_ele,C_ele,s);               % C(s)
    tmp02=interp1(s_ele,R51_ele,s);             % R51(s)
    tmp03=interp1(s_ele,R52_ele,s);             % R52(s)
    tmp04=interp1(s_ele,R53_ele,s);             % R53(s)
    tmp05=interp1(s_ele,R54_ele,s);             % R54(s)
    tmp06=interp1(s_ele,R56_ele,s);             % R56(s)
    
    if (iIBS==0)
        tmp07=(tmp01.^2)*(k_wave^2).*(emitx)/(2*betax0);
    else
        tmp07=(tmp01.^2)*(k_wave^2).*(emit_gx_IBS)/(2*betax0);
    end
    
    %tmp08=(betax0^2)*(tmp02.^2)+(tmp03.^2);
    tmp08=(betax0^2)*((tmp02-tmp03*alphax0/betax0).^2)+(tmp03.^2);
    
    if (iIBS==0)
        tmp09=(tmp01.^2)*(k_wave^2).*(emity)/(2*betay0);
    else
        tmp09=(tmp01.^2)*(k_wave^2).*(emit_gy_IBS)/(2*betay0);
    end
    
    %tmp10=(betay0^2)*(tmp04.^2)+(tmp05.^2);
    tmp10=(betay0^2)*((tmp04-tmp05*alphay0/betay0).^2)+(tmp05.^2);
    
    if (iIBS==0)
        tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    else
        tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta_IBS.^2).*(tmp06.^2)/2;
    end
    
    LD=exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    
    func=n_1k0*LD;
    %---------------------------------------------------------------------%
    
    
    
    
    
    