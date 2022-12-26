%function func=call_dipole

global dipole_s N_D;

load('lattice_dipole_info.o');

dipole_s=lattice_dipole_info*100; % unit: cm
[r,~]=size(dipole_s);

N_D=r/2; % determine the number of dipoles

if (rem(N_D,1)~=0)
    fprintf('warning: not an integer number of dipoles from input lattice... \n');
    %pause;
end

