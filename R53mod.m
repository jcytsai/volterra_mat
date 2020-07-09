function func=R53mod(a,b)
    global s_ele;
    global R53_ele;
    global C_ele;
    tmp01=interp1(s_ele,C_ele,a);       % C(s)
    tmp02=interp1(s_ele,R53_ele,a);     % R53(s)
    tmp03=interp1(s_ele,C_ele,b);       % C(s')
    tmp04=interp1(s_ele,R53_ele,b);     % R53(s')
    func=tmp01.*tmp02-tmp03.*tmp04;