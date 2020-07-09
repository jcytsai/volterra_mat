function func=R51mod(a,b)
    global s_ele;
    global R51_ele;
    global C_ele;
    tmp01=interp1(s_ele,C_ele,a);       % C(s)
    tmp02=interp1(s_ele,R51_ele,a);     % R51(s)
    tmp03=interp1(s_ele,C_ele,b);       % C(s')
    tmp04=interp1(s_ele,R51_ele,b);     % R51(s')
    func=tmp01.*tmp02-tmp03.*tmp04;