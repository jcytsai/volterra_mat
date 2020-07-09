function func=csr1d(k,s,index_tr)

A=3^(-1/3)*gamma(2/3)*(1i*sqrt(3)-1);
tmp01=auxr(s);

%if (index_tr==0)
    if (tmp01 > 1e10)
        func=0.0;
    else
        %if (s<4730)
        %    func=0.0;
        %else
        func=-1i*k^(1/3)*A*abs(tmp01)^(-2/3);
        %end
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
            func=-1i*k^(1/3)*A*abs(tmp01)^(-2/3);
        end
    else
        func=0.0;
    end
end
%}
