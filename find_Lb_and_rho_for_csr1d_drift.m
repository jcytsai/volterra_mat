function [func1,func2]=find_Lb_and_rho_for_csr1d_drift(t)

global dipole_s N_D;

upstream_dipole_Lb=0.0;
upstream_dipole_rho=1e50;

for m=1:1:N_D
   if (t>dipole_s(2*m,1))
       upstream_dipole_Lb(m)=abs(dipole_s(2*m,1)-dipole_s(2*m-1,1));
       if (abs(dipole_s(2*m-1,2)) < abs(dipole_s(2*m-1,3))) % horizontal bend
           upstream_dipole_rho(m)=dipole_s(2*m-1,2);
       else                                                 % vertical bend
           upstream_dipole_rho(m)=dipole_s(2*m-1,3);
       end
   end
end

[r,c]=size(upstream_dipole_Lb);
if (c==0)
    upstream_dipole_Lb(1)=0.0; % i.e. no dipole upstream the beamline
    upstream_dipole_rho(1)=1e50;
end

func1=upstream_dipole_Lb';
func2=upstream_dipole_rho';