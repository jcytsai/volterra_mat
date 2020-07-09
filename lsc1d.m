function func=lsc1d(k,Sx,Sy,s,model)
    global s_ele egamma_vec nb round_pipe_radius;
    
    % the following constant quantities are in MKS units
    electron_mass=9.10938291e-31;
    electron_charge=1.6e-19;
    epsilon_0=8.854187817e-12;
    
    c_speed=2.99792458e10;
    
    
    egamma=interp1(s_ele,egamma_vec,s);
    rb=1.747/2*(Sx+Sy);    % unit: cm
    %rb=0.05;
    xi=k*rb/egamma;
    
    %omega_pe=sqrt(electron_charge^2*nb/(egamma^3*epsilon_0*electron_mass));
    %plasma_osc=cos(omega_pe/c_speed*s);
    
    if (model==1)
    % model A, on-axis (transverse uniform, circular cross section)
        tmp01=4*1i/rb/egamma;
        tmp02=(1-xi*besselk(1,xi))/xi;
        func=tmp01*tmp02;
    end
    
    if (model==2)
    % model B, average (averaged over transverse beam distribution, for Gaussian or parabolic)
        tmp01=4*1i/rb/egamma;
        tmp02=(1-2*besseli(1,xi).*besselk(1,xi))./xi;
        func=tmp01*tmp02;
    end
    
    if (model==3)
    % model C, transverse axisymmetric Gaussian model
        S_x_y_ave=0.5*(Sx+Sy);
        xi_a=k*S_x_y_ave/egamma;
        tmp01=xi_a/(egamma*S_x_y_ave);
        tmp02=0.5*xi_a^2;
        tmp03=exp(tmp02);
        tmp04=-expint(tmp02);    % Note: Matlab definition is different from that in PRST-AB 11, 034401 (2008)       
        func=-1i*tmp01*tmp03*tmp04;
    end
    
    if (model==4)
    % model D, on-axis (transverse uniform, circular cross section) with
    % perfectly conducting beam pipe rp
        xii=xi*round_pipe_radius/rb;
        tmp01=4*1i/rb/egamma;
        tmp02=(1-xi*(besselk(1,xi)+besselk(0,xii)*besseli(1,xi)/besseli(0,xii)))/xi;
        func=tmp01*tmp02;
    end
    
    %fprintf('LSC impedance Z=%f+%fi at s=%f...\n',real(func),imag(func),s/100);
    