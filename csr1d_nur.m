function func=csr1d_nur(k,s,egamma,index_tr,index_switch)

func1=@(x) airy(0,x);                                 % integrand of Eq.(8)
func2=@(t) airy(2,t);                                 % first integrand of Eq.(12), Bi(t)
func3=@(t) airy(0,t);                                 % second integrand of Eq.(12), Ai(t)
func4=@(x) besselk(5/3,x);
func5=@(t,x) t.*cos(t.^3/3+x.*t); 

tmp01=abs(auxr(s));

%if (index_tr==0)
if (tmp01 > 1e10)
    func=0.0;
else
    
    tmp02=k^(1/3)*abs(tmp01)^(-2/3);
    tmp03=(k*tmp01)^(2/3)/egamma^2;
    tmp04=quadgk(func1,0,tmp03)-1/3;
    tmp05=2*k*tmp01/(3*egamma^3);
    tmp06=quadgk(func2,0,tmp05);
    tmp07=quadgk(func3,0,tmp05);
    
    Re_Z=tmp02*(-2*pi)*airy(1,tmp03)+k*pi/egamma^2*tmp04;
    Im_Z=tmp02*(2*pi/3*airy(3,tmp05)+2*pi*airy(1,tmp05)*tmp06-2*pi*airy(3,tmp05)*tmp07);

    %{
    tmp08=k/(egamma^2*sqrt(3));
    tmp09=k*tmp01/(3*egamma^3);
    tmp10=quadgk(func4,2*tmp09,Inf);
    Re_Z=tmp08*tmp10;
    tmp11=2*k^3/tmp01^(2/3);
    func6=@(t) func5(t,tmp03);
    tmp12=quadgk(func6,0,Inf);
    Im_Z=tmp11*tmp12;
    %}
    
    func=Re_Z+1i*Im_Z;
end
%{
else
    local_s=find_local_s(s);
    lambda_tmp=2*pi/k;
    formation_length=(24*lambda_tmp*tmp01^2)^(1/3);
    if (local_s > formation_length) 
        if (tmp01 > 1e10)
            func=0.0;
        else
    
            tmp02=k^(1/3)*abs(tmp01)^(-2/3);
            tmp03=(k*tmp01)^(2/3)/egamma^2;
            tmp04=quadgk(func1,0,tmp03)-1/3;
            tmp05=2*k*tmp01/(3*egamma^3);
            tmp06=quadgk(func2,0,tmp05);
            tmp07=quadgk(func3,0,tmp05);
    
            Re_Z=tmp02*(-2*pi)*airy(1,tmp03)+k*pi/egamma^2*tmp04;
            Im_Z=tmp02*(2*pi/3*airy(3,tmp05)+2*pi*airy(1,tmp05)*tmp06-2*pi*airy(3,tmp05)*tmp07);
    
            func=Re_Z+1i*Im_Z;
        end
    else
        func=0.0;
    end
end
%}

