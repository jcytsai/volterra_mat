function func=find_local_s(t)

global dipole_s N_D;

func=0.0;
for m=1:1:N_D
   if (t>=dipole_s(2*m-1,1) && t<=dipole_s(2*m,1))
       func=t-dipole_s(2*m-1,1);
   end
end
