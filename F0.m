function func=F0(bbeta)

%bbeta2=bbeta.^2;
%tmp01=airy(1,bbeta2).*(airy(0,bbeta2)-1i*airy(3,bbeta2));
%tmp02=bbeta2.*airy(0,bbeta2).*(airy(0,bbeta2)-1i*airy(2,bbeta2));
%func=tmp01+tmp02;

if (bbeta>10)
    tmp01=bbeta/(2*pi)*exp(-4/3*bbeta^3);
    tmp02=1+1/(24*bbeta^3)+1/(1152*bbeta^6);
    tmp03=-3/(16*pi*bbeta^5);
    tmp04=1+105/32/bbeta^6;
    Re_F0=tmp01*tmp02;
    Im_F0=tmp03*tmp04;
    func=Re_F0+1i*Im_F0;
else
    bbeta2=bbeta.^2;
    tmp01=airy(3,bbeta2); %if (tmp01 > 1e300); tmp01=1e300; end
    tmp02=airy(2,bbeta2); %if (tmp02 > 1e300); tmp02=1e300; end
    tmp03=airy(1,bbeta2).*airy(1,bbeta2)-1i*airy(1,bbeta2).*tmp01;
    tmp04=bbeta2.*airy(0,bbeta2).*airy(0,bbeta2)-1i*bbeta2.*airy(0,bbeta2).*tmp02;
    func=sum(tmp03+tmp04);
end

