function func=R52mod(a,b)
    global s_ele;
    global R52_ele;
    global C_ele;
    tmp01=interp1(s_ele,C_ele,a);       % C(s)
    tmp02=interp1(s_ele,R52_ele,a);     % R52(s)
    tmp03=interp1(s_ele,C_ele,b);       % C(s')
    tmp04=interp1(s_ele,R52_ele,b);     % R52(s')
    func=tmp01.*tmp02-tmp03.*tmp04;