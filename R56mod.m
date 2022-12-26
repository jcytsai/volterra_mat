function func=R56mod(a,b)
    global s_ele;
    global R56_ele;
    global C_ele;
    tmp01=interp1(s_ele,C_ele,a);          % C(s)
    tmp02=interp1(s_ele,R56_ele,a);        % R56(s)
    tmp03=interp1(s_ele,C_ele,b);          % C(s')
    tmp04=interp1(s_ele,R56_ele,b);        % R56(s')
    func=tmp01.*tmp02-tmp03.*tmp04;