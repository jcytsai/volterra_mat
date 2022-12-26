function func=Piwinski_g_func(w)
% w: an array of the size same as s_array
% return an array

z_array=(w.^2+1)./(2*w);
func_01=@(x,z) 1./sqrt(z+sqrt(z^2-1)*cos(x))/pi;  % z: scalar parameter

for m=1:1:length(w)
    if (z_array(m)==1)
        P_func_01(m)=1;
    else
        P_func_01(m)=integral(@(x)func_01(x,z_array(m)),0,pi);
    end
    %if (isnan(P_func_01(m)))
    %    pause;
    %end
end

for m=1:1:length(w)
    if (z_array(m)==1)
        P_func_02(m)=0;
    else
        tmp01=1./sqrt((z_array(m))^2-1);
        if (abs(tmp01)>1e4)
            fprintf('tmp01=%f too large, check Piwinski_g_func.m...\n',abs(tmp01));
        end
        if (m==1)
            Y=[1 P_func_01(1)];
            X=linspace(1,z_array(m),2);
        else
            Y=P_func_01(1:m);
            X=linspace(1,z_array(m),m);
        end
        P_func_02(m)=tmp01*trapz(X,Y);
    end
end

%P_func_02(1)=interp1(w(2:end),P_func_02(2:end),w(1),'pchip','extrap');
func=sqrt(pi./w).*(P_func_01+1.5*sign(w-1).*P_func_02);

% 2020.02.25, this subroutine was benchmarked against fermilab-pub-89-224,
% Fig. 3 and 4 with the following commands:
% y=Piwinski_g_func(b);
% yy=Piwinski_g_func(1./b);
% plot(b,y);
% plot(b,y+yy./b);