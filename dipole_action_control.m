function func=dipole_action_control(t,col,iflag)

global dipole_control_seq dipole_s N_D

func=iflag;
for m=1:1:N_D
   if (t<=dipole_s(2*m,1) && t>=dipole_s(2*m-1,1))
       func=dipole_control_seq(m,col);
       if (func~=iflag && iflag==1)
           fprintf('CSR off in %d-th dipole...\n',m);
       end
   end
end