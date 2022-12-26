function func=g0kyz_mat(s)

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
    global alphax0;
    global alphay0;
    global betax0;
    global betay0;
    global gammax0;
    global gammay0;
    global sigma_delta;
    global n_1k0;
    global e_1k0;
    global ax_1k0 axp_1k0 ay_1k0 ayp_1k0;
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
    %tmp08=(betax0^2)*(tmp02.^2)+(tmp03.^2);
    tmp08=(betax0^2)*((tmp02-tmp03*alphax0/betax0).^2)+(tmp03.^2);
    tmp09=(tmp01.^2)*(k_wave^2).*(emity)/(2*betay0);
    %tmp10=(betay0^2)*(tmp04.^2)+(tmp05.^2);
    tmp10=(betay0^2)*((tmp04-tmp05*alphay0/betay0).^2)+(tmp05.^2);
    
    %if (find_RF_para==1)
        %tmp11=(tmp01.^2).*(k_wave^2).*((sigma_delta*tmp_egamma_vec).^2).*(tmp06.^2)/2;
        tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    %else
    %    tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    %end
    
    tmp12=2*alphay0*gammay0*tmp05;
    tmp13=(betay0*gammay0+alphay0^2)*tmp04;
    
    LD=exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    
    func=1i*k_wave*tmp01.*(tmp12-tmp13)*ay_1k0.*LD;
    