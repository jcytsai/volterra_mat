function func=R56transport(b,a)
% this version modifies R56transport function such that it can adopt the
% second argument a as a vector while keeping the first argument b as a
% scalar
% output of this subroutine: a vector func
%
% b=s'(previous time), a=s(present time)
%
    global s_ele R51_ele R52_ele R53_ele R54_ele R56_ele R55_ele;
    
    %tmp01=interp1(s_ele,R56_ele,a); % vector
    %tmp02=interp1(s_ele,R56_ele,b);
    tmp01=interp1(s_ele,R56_ele,a).*interp1(s_ele,R55_ele,b); % vector
    tmp02=interp1(s_ele,R56_ele,b).*interp1(s_ele,R55_ele,a); % vector
    tmp03=interp1(s_ele,R51_ele,b).*interp1(s_ele,R52_ele,a); % vector
    tmp04=interp1(s_ele,R51_ele,a).*interp1(s_ele,R52_ele,b); % vector
    tmp05=interp1(s_ele,R53_ele,b).*interp1(s_ele,R54_ele,a); % vector
    tmp06=interp1(s_ele,R53_ele,a).*interp1(s_ele,R54_ele,b); % vector
    
    func=(tmp01-tmp02+tmp03-tmp04+tmp05-tmp06); % vector
    
    
    