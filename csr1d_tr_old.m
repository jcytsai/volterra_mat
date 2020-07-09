function func=csr1d_tr(k,s)

local_s=find_local_s(s);            % unit: cm
if (local_s < 0.1); local_s=0.1; end

if (local_s==0)
    func=0.0;
else
    
z_L=local_s^3/(24*auxr(s)^2);

%if ((z_L<1e-7) || (local_s<1)); z_L=1e-7; local_s=1; end

method=2;
% D. Zhou's version
if (method==1)
K=-2/(3^(1/3)*abs(auxr(s))^(2/3));
gamma_DZ=(1-gammaincc(1i*k*z_L,-1/3))*gamma(-1/3);
tmp01=K*exp(-1i*4*k*z_L)/z_L^(1/3);
tmp02=-1/3*K*(1i*k)^(1/3)*gamma_DZ;
tmp03=tmp01+tmp02;
else
mu=k*z_L;
tmp01=-4/local_s*exp(-1i*4*mu);
tmp02=4/(3*local_s)*(1i*mu)^(1/3)*mfun('GAMMA',-1/3,1i*mu);
tmp03=tmp01+tmp02;
end

%{
% C. Mitchell & J. Qiang's version
mu=k*z_L;
tmp01=exp(-1i*mu)-exp(-4*1i*mu);
%tmp01=-exp(-4*1i*mu);
%gamma_MQ=(1-gammaincc(1i*mu,2/3))*gamma(2/3);
tmp02=-(1i*mu)^(1/3)*mfun('GAMMA',2/3,1i*mu);
tmp03=4/local_s*(tmp01+tmp02);
%}

%if (abs(tmp03)>10) tmp03=0.0; fprintf('warning...\n'); end

func=tmp03;
end