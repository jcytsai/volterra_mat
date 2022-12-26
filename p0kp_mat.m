function func=p0kp_mat(s)

    global s_ele;
    global R51_ele;
    global R52_ele;
    global R53_ele;
    global R54_ele;
    global R56_ele;
    global C_ele;
    global k_wave;
    global emitx;
    global emity;
    global betax0;
    global betay0;
    global sigma_delta;
    global n_1k0;
    global e_1k0;
    %global egamma_vec emit_norm_x emit_norm_y find_RF_para;
    
    format long
    
    %tmp_egamma_vec=interp1(s_ele,egamma_vec,s);
    %emitx=emit_norm_x./tmp_egamma_vec;
    %emity=emit_norm_y./tmp_egamma_vec;
    
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
    tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    
    %func=-1i*k_wave*tmp01.*tmp06*sigma_delta^2*n_1k0.*exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    func=(1-k_wave^2*tmp01.^2.*tmp06.^2*sigma_delta^2)*e_1k0.*exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);