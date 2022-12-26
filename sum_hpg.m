function func=sum_hpg(N_max,arg_03,sum_order)

tmp_sum=zeros(size(arg_03));

if (sum_order==1.5)
    if (N_max==0)
        tmp02=gamma(3/2);
        tmp03=cos(sqrt(abs(-arg_03)*3));     % vector, much faster
        tmp_sum=tmp02*tmp03;
    else
    for n=0:1:N_max
        if (n==0)
            tmp01=1;
        else
            tmp01=(-2)^n*prod(1:2:2*n-1)/prod(1:1:2*n+1);
        end
        tmp02=gamma(n+3/2);
        %tmp03=hypergeomq(n+3/2,1/2,arg_03);       % vector, fast
        %tmp03=hypergeom(n+3/2,1/2,arg_03);        % vector, slow
        tmp03=cos(sqrt(abs(-arg_03)*(2*n+3)));     % vector, much faster
        tmp_sum=tmp_sum+tmp01*tmp02*tmp03;
    end
    end
    func=tmp_sum;
end

if (sum_order==2.5)
    if (N_max==0)
        tmp02=gamma(5/2);
        tmp03=sinc(sqrt(abs(-arg_03)*5)/pi); % vector, much faster
        tmp_sum=tmp02*tmp03;
    else
    for n=0:1:N_max
        if (n==0)
            tmp01=1;
        else
            tmp01=(-2)^n*prod(1:2:2*n-1)/prod(1:1:2*n+1);
        end
        tmp02=gamma(n+5/2);
        %tmp03=hypergeomq(n+5/2,3/2,arg_03);           % vector, fast
        %tmp03=hypergeom(n+5/2,3/2,arg_03);            % vector, slow
        tmp03=sinc(sqrt(abs(-arg_03))*sqrt(2*n+5)/pi); % vector, much faster
        tmp_sum=tmp_sum+tmp01*tmp02*tmp03;
    end
    end
    func=tmp_sum;
end
