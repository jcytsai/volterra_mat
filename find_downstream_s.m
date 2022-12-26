function func=find_downstream_s(t)

global dipole_s N_D;

tmp_s=0.0;
for m=1:1:N_D
   if (t>dipole_s(2*m,1))
       tmp_s(m)=t-dipole_s(2*m,1);
   end
end

[r,c]=size(tmp_s);
if (c==0); tmp_s(1)=0.0; end

func=tmp_s';

