function func=ay0kp_mat(s)

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
    global R11_ele R12_ele R13_ele R14_ele R21_ele R22_ele R23_ele R24_ele R31_ele R32_ele R33_ele R34_ele R41_ele R42_ele R43_ele R44_ele R16_ele R26_ele R36_ele R46_ele;
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
    tmp11=(tmp01.^2).*(k_wave^2).*(sigma_delta^2).*(tmp06.^2)/2;
    tmp12=tmp03*alphax0-tmp02*betax0;           % R52(s)*alphax0-R51(s)*betax0
    tmp13=tmp03*gammax0-tmp02*alphax0;          % R52(s)*gammax0-R51(s)*alphax0
    tmp14=tmp05*alphay0-tmp04*betay0;           % R54(s)*alphay0-R53(s)*betay0
    tmp15=tmp05*gammay0-tmp04*alphay0;          % R54(s)*gammay0-R53(s)*alphay0
    %tmp16=interp1(s_ele,R11_ele,s);             % R11(s);
    %tmp17=interp1(s_ele,R12_ele,s);             % R12(s);
    %tmp18=interp1(s_ele,R13_ele,s);             % R13(s);
    %tmp19=interp1(s_ele,R14_ele,s);             % R14(s);
    %tmp20=interp1(s_ele,R21_ele,s);             % R21(s);
    %tmp21=interp1(s_ele,R22_ele,s);             % R22(s);
    %tmp22=interp1(s_ele,R23_ele,s);             % R23(s);
    %tmp23=interp1(s_ele,R24_ele,s);             % R24(s);
    tmp24=interp1(s_ele,R31_ele,s);             % R31(s);
    tmp25=interp1(s_ele,R32_ele,s);             % R32(s);
    tmp26=interp1(s_ele,R33_ele,s);             % R33(s);
    tmp27=interp1(s_ele,R34_ele,s);             % R34(s);
    %tmp28=interp1(s_ele,R41_ele,s);             % R41(s);
    %tmp29=interp1(s_ele,R42_ele,s);             % R42(s);
    %tmp30=interp1(s_ele,R43_ele,s);             % R43(s);
    %tmp31=interp1(s_ele,R44_ele,s);             % R44(s);
    %tmp32=interp1(s_ele,R16_ele,s);             % R16(s);
    %tmp33=interp1(s_ele,R26_ele,s);             % R26(s);
    tmp34=interp1(s_ele,R36_ele,s);             % R36(s);
    %tmp35=interp1(s_ele,R46_ele,s);             % R46(s);
    
    LD=exp(-(tmp07).*(tmp08)-(tmp09).*(tmp10)-tmp11);
    
    %func01=1i*(k_wave*tmp01).*tmp12*emitx*n_1k0.*LD;           % x0,d
    func02=(k_wave*tmp01).^2*emitx.*tmp06.*tmp12*e_1k0.*LD;    % x0,e
    %func03=1i*(k_wave*tmp01).*tmp13*emitx*n_1k0.*LD;           % x0',d
    func04=(k_wave*tmp01).^2*emitx.*tmp06.*tmp13*e_1k0.*LD;    % x0',e
    
    %func05=1i*(k_wave*tmp01).*tmp14*emity*n_1k0.*LD;           % y0,d
    func06=(k_wave*tmp01).^2*emity.*tmp06.*tmp14*e_1k0.*LD;    % y0,e
    %func07=1i*(k_wave*tmp01).*tmp15*emity*n_1k0.*LD;           % y0',d
    func08=(k_wave*tmp01).^2*emity.*tmp06.*tmp15*e_1k0.*LD;    % y0',e
    
    %func09=-1i*(k_wave*tmp01).*tmp06*sigma_delta^2*n_1k0.*LD;  % delta0,d
    func10=(1-k_wave^2*tmp01.^2.*tmp06.^2*sigma_delta^2)*e_1k0.*LD; % delta0,e
    
    func=tmp24.*func02+tmp25.*func04+tmp26.*func06+tmp27.*func08+tmp34.*func10;
    