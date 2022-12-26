function func=auxr(t)
% this subroutine outputs rho_x and rho_y with sign and magnitude in both
% horizontal and vertical planes; outputs INFINITY if neither of both cases

global dipole_s N_D;

func=1e50;

for m=1:1:N_D
%if (N_D >= m)
    if (t>=dipole_s(2*m-1,1) && t<=dipole_s(2*m,1))
        if (abs(dipole_s(2*m-1,2)) < abs(dipole_s(2*m-1,3))) % horizontal bend
            func=dipole_s(2*m-1,2);
        else                                                 % vertical bend
            func=dipole_s(2*m-1,3);
        end
    end
%end
end
