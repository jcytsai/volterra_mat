function func=find_Lb(t)

global dipole_s N_D;

for m=1:1:N_D
   if (t>=dipole_s(2*m-1,1) && t<=dipole_s(2*m,1))
       func=abs(dipole_s(2*m,1)-dipole_s(2*m-1,1));
   end
end
