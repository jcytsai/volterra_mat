function func=R54transport(b,a)
% this version modifies R56transport function such that it can adopt the
% second argument a as a vector while keeping the first argument b as a
% scalar 
% output of this subroutine: a vector func
%
% b=s' or tau(previous time,scalar), a=s(present time,vector)
%
    global s_ele R53_ele R54_ele R33_ele R34_ele;
    
    tmp01=interp1(s_ele,R33_ele,b);                           % scalar
    tmp02=interp1(s_ele,R34_ele,b);                           % scalar
    tmp03=interp1(s_ele,R54_ele,a)-interp1(s_ele,R54_ele,b);  % vector
    tmp04=interp1(s_ele,R53_ele,a)-interp1(s_ele,R53_ele,b);  % vector
    
    func=tmp03*tmp01-tmp04*tmp02; % vector