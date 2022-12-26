function func=C_vec(s)

    global s_ele R51_ele R52_ele R53_ele R54_ele R56_ele C_ele k_wave betax0 betay0 sigma_delta emit_norm_x emit_norm_y egamma_vec;
    
    format long
    
    emitx=interp1(s_ele,emit_norm_x./egamma_vec,s);
    emity=interp1(s_ele,emit_norm_y./egamma_vec,s);    
    
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
    
    func=1i*k_wave*tmp01.*tmp06.*exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    