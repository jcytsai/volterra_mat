function func=R26transport(b,a)
%
% b=s'(previous time, scalar), a=s(present time, vector)
%
    global s_ele R11_ele R12_ele R13_ele R14_ele R21_ele R22_ele R23_ele R24_ele R31_ele R32_ele R33_ele R34_ele R41_ele R42_ele R43_ele R44_ele R16_ele R26_ele R36_ele R46_ele R51_ele R52_ele R53_ele R54_ele R55_ele R56_ele;
    
    tmp01=interp1(s_ele,R22_ele,a).*interp1(s_ele,R51_ele,b); % vector
    tmp02=interp1(s_ele,R52_ele,b).*interp1(s_ele,R21_ele,a); % vector
    tmp03=interp1(s_ele,R26_ele,a);                           % vector
    
    func=tmp03+tmp01-tmp02;                                   % vector