% This program intends to combine different sections of lattice and to
% append the transport functions such that they are continuous along the
% whole lattice for subsequent volterra_mat simulation.

clear all
format long

%---------------------- required input parameters ------------------------%
rep_s_start=6.518;      % unit: m
chirp01=-21.4;          % unit: m^(-1)
chirp02=39.83;
%-------------------------------------------------------------------------%

% load transport functions from assigned sections of lattices
% foramt [s_ele R16_ele R36_ele R51_ele R52_ele R53_ele R54_ele R55_ele
% R56_ele egamma_vec C_factor] in MKS unit
filename01='lattice_transport_functions_ELEGANT_corrected.o';
filename02='lattice_transport_functions_rep.o';

delimiterIn=' '; headerlinesIn=0;

transport01=importdata(filename01,delimiterIn,headerlinesIn);
transport02=importdata(filename02,delimiterIn,headerlinesIn);


% append transport functions
rep_lattice_length=transport02(end,1)-transport02(1,1);

rep_s_end=transport02(end,1)+rep_s_start;

[r01 c01]=size(transport01);
[r02 c02]=size(transport02);

% use m to determine starting row index
for m=1:1:r01
    if (transport01(m,1)>=rep_s_start)
        break
    end
end


% use n to determine ending row index
for n=m:1:r01
    if (transport01(n,1)>=rep_s_end)
        break
    end
end

% merge and generate/output a new lattice_transport_functions_combined.o file
new_egamma_vec=[ones(m-1,1)*transport02(1,end-1);transport02(:,end-1);ones(n:end,1)*transport02(end,end-1)];
[rnew cnew]=size(new_egamma_vec);


R55_01=transport01(1:m-1,8);
R56_01=transport01(1:m-1,9);

if (R55_01(1)~=1)
    R55_01=R55_01-R55_01(1);
end
if (R56_01(1)~=0)
    R56_01=R56_01-R56_01(1);
end

C_factor_01=1./(R55_01-chirp01*R56_01);

R55_02=transport01(n:end,8);
R56_02=transport01(n:end,9);

if (R56_02(1)~=0)
    R56_02=R56_02-R56_02(1);
    R56_02=-R56_02;
end

C_factor_02=1./(1-chirp02*R56_02);

new_C_factor=[C_factor_01;transport02(:,end);C_factor_02];

new_s_ele=[transport01(1:m-1,1);transport02(:,1)+rep_s_start;transport01(m+r02:end,1)];
new_R_mat=[transport01(1:m-1,2:9);(transport02(:,2:9)-ones(r02,1)*transport01(m-1,2:9));transport01(m+r02:end,2:9)];
new_transport=[new_s_ele new_R_mat new_egamma_vec new_C_factor];

pause;


