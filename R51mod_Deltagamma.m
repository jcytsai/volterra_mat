function func=R51mod_Deltagamma(a,b)
    global s_ele;
    global R51_ele;
    global C_ele;
    global egamma_vec;
    
    egamma_a=interp1(s_ele,egamma_vec,a);            % vector
    egamma_b=interp1(s_ele,egamma_vec,b);            % scalar
    
    tmp01=interp1(s_ele,C_ele,a);                    % C(s)
    tmp02=interp1(s_ele,R51_ele,a).*egamma_a;        % R51(s)
    tmp03=interp1(s_ele,C_ele,b);                    % C(s')
    tmp04=interp1(s_ele,R51_ele,b)*egamma_b;         % R51(s')
    func=tmp01.*tmp02-tmp03.*tmp04;