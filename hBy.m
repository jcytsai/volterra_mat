function [func,rhoy_s]=hBy(s_array)
% this subroutine returns 1 if s is located within a dipole and 0 elsewhere
% s is an input array in unit of cm; output is an array of the same size

global dipole_s N_D;

[~,c]=size(s_array);

func=zeros(1,c);
rhoy_s=1e50*ones(1,c);

for p=1:1:c
    s=s_array(p);
    for m=1:1:N_D
        if (s>=dipole_s(2*m-1,1) && s<=dipole_s(2*m,1))
            if (abs(dipole_s(2*m-1,2)) > abs(dipole_s(2*m-1,3))) % vertical bend
                func(p)=1;
                rhoy_s(p)=abs(dipole_s(2*m-1,2));
            end
        end
    end
end

