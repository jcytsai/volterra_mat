%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subroutine reads several relevant input parameters for volterra_mat%
% from ELEGANT input/output files and from user inputs.                   %
% Estimation of bunch current assumes Gaussian bunch distribution, i.e.   %
% I_b=q_charge/2.35*sigma_t                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
clearvars -EXCEPT shortest;
format long

global energy I_b C_factor emit_norm_x emit_norm_y betax0 betay0 alphax0 alphay0 sigma_delta chirp start_pos end_pos lambda_start01 lambda_end01 scan_num01 mesh_num lambda_start02 lambda_end02 scan_num02 iplot_lattice iplot_gain_func iplot_gain_spec isave;
global iiterative ilast iCSR_ss iCSR_tr iCSR_drift iLSC ilinac iGaussian iEnergy_mod_calc iEnergy_gain_calc iplot_energy_mod surfplot_gain_spec_func ssCSR_model LSC_model issCSRpp iquilt_plot idm_analysis iplot_I_b iRF_ele find_TWLA full_pipe_height;
global Derbenev iTransverse_gain_calc first_harmonic_notification round_pipe_radius;

skip_ask=1;

%prompt='use GUI [1] or direct run [2]...';
%run_mode=input(prompt);

disp('Read input parameters from ELEGANT input/output files...');

filename01='rootname.ele';
filename02='rootname.lte';
filename03='lattice_transport_functions.o';

fid =fopen(filename01,'rt');	%open the file
fid2=fopen(filename02,'rt');	%open the file

%-------------------------------------------------------------------------%
line=' ';
% load beam energy in GeV
while (isempty(strfind(line,'p_central_mev')) && isempty(strfind(line,'p_central')))
  line=fgetl(fid);
end
energy_MeV=sscanf(line,' p_central_mev = %f');
energy_bg=sscanf(line,' p_central = %f');

if (isempty(energy_MeV)==0)
    energy=energy_MeV/1e3;
else
    energy=(0.51099891/1e3)*sqrt(energy_bg^2+energy_bg*sqrt(energy_bg^2-4))/sqrt(2);
end

if ((isempty(energy_bg)==1) && (isempty(energy_MeV)==1))
    fprintf('neither p_central nor p_central_mev given... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
line=' ';
% load beam emittances
while (isempty(strfind(line,'emit_nx')) && isempty(strfind(line,'emit_x')))
  line=fgetl(fid);
end
emit_norm_x=sscanf(line,' emit_nx = %f')*1e6;
emit_geom_x=sscanf(line,' emit_x = %f')*1e6;

if (isempty(emit_norm_x)==0)
    emit_norm_x=emit_norm_x;
else
    if (isempty(energy_MeV)==0)
        eegamma=energy_MeV/(0.51099891);
        emit_norm_x=emit_geom_x*(eegamma/sqrt(1-1/eegamma^2));
    else
        emit_norm_x=emit_geom_x*energy_bg;
    end
end

if ((isempty(emit_norm_x)==1) && (isempty(emit_geom_x)==1))
    fprintf('neither emit_geom_x nor emit_norm_x given... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
line=' ';
while (isempty(strfind(line,'emit_ny')) && isempty(strfind(line,'emit_y')))
  line=fgetl(fid);
end
emit_norm_y=sscanf(line,' emit_ny = %f')*1e6;
emit_geom_y=sscanf(line,' emit_y = %f')*1e6;

if (isempty(emit_norm_y)==0)
    emit_norm_y=emit_norm_y;
else
    if (isempty(energy_MeV)==0)
        eegamma=energy_MeV/(0.51099891);
        emit_norm_y=emit_geom_y*(eegamma/sqrt(1-1/eegamma^2));
    else
        emit_norm_y=emit_geom_y*energy_bg;
    end
end

if ((isempty(emit_norm_y)==1) && (isempty(emit_geom_y)==1))
    fprintf('neither emit_geom_y nor emit_norm_y given... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load initial beta functions
line=' ';
while (isempty(strfind(line,'beta_x')))
  line=fgetl(fid);
end
betax0=sscanf(line,' beta_x = %f');

if ((isempty(betax0)==1))
    betax0=1e-8;
    fprintf('betax0 not given, set default 0... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file



line=' ';
while (isempty(strfind(line,'beta_y')))
  line=fgetl(fid);
end
betay0=sscanf(line,' beta_y = %f');

if ((isempty(betay0)==1))
    betay0=1e-8;
    fprintf('betax0 not given, set default 0... \n');
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load initial alpha functions
while (isempty(strfind(line,'alpha_x')))
  line=fgetl(fid);
end
alphax0=sscanf(line,' alpha_x = %f');

if ((isempty(alphax0)==1))
    alphax0=1e-8;
    fprintf('alphax0 not given, set default 0... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file



while (isempty(strfind(line,'alpha_y')))
  line=fgetl(fid);
end
alphay0=sscanf(line,' alpha_y = %f');

if ((isempty(alphay0)==1))
    alphay0=1e-8;
    fprintf('alphay0 not given, set default 0... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load initial uncorrelated energy spread
while (isempty(strfind(line,'sigma_dp')))
  line=fgetl(fid);
end
sigma_delta=sscanf(line,' sigma_dp = %f');

if ((isempty(sigma_delta)==1))
    sigma_delta=1e-30;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load bunch length (assume Gaussian)
while (isempty(strfind(line,'sigma_s')))
  line=fgetl(fid);
end
sigma_s=sscanf(line,' sigma_s = %f');
sigma_t=sigma_s/(2.99792458e8);
fwhm_t=2.35*sigma_t;  % unit: sec

if ((isempty(sigma_s)==1))
    fprintf('sigma_s not given... \n'); %pause;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load chirp parameter
while (isempty(strfind(line,'momentum_chirp')))
  line=fgetl(fid);
  if (line==-1)
      break;
  end
end
if (line~=-1)
    chirp=sscanf(line,' momentum_chirp = %f');
else
    chirp=0;
end
fclose(fid);
fid  = fopen(filename01,'rt');	%open the file
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load bunch charge
line=' '; loop_num=0;
while (isempty(strfind(line,'q: charge,total')) && isempty(strfind(line,'q: charge, total')) && isempty(strfind(line,'Q: CHARGE,TOTAL')) && isempty(strfind(line,'Q: CHARGE, TOTAL')) && isempty(strfind(line,'Q: charge,total')) && isempty(strfind(line,'Q: charge, total')) && loop_num <= 500)
  line=fgetl(fid2); loop_num=loop_num+1;
end

if (line~=-1)
    if (isempty(strfind(line,'q: charge,total'))==0)
        q_charge=sscanf(line,'q: charge,total=%f'); % unit: C
    elseif (isempty(strfind(line,'q: charge, total'))==0)
        q_charge=sscanf(line,'q: charge, total=%f'); % unit: C
    elseif (isempty(strfind(line,'Q: CHARGE,TOTAL'))==0)
        q_charge=sscanf(line,'Q: CHARGE,TOTAL=%f'); % unit: C
    elseif (isempty(strfind(line,'Q: CHARGE, TOTAL'))==0)
        q_charge=sscanf(line,'Q: CHARGE, TOTAL=%f'); % unit: C
    elseif (isempty(strfind(line,'Q: charge,total'))==0)
        q_charge=sscanf(line,'Q: charge,total=%f'); % unit: C
    elseif (isempty(strfind(line,'Q: charge, total'))==0)
        q_charge=sscanf(line,'Q: charge, total=%f'); % unit: C
    end
end

if (line==-1)
    fprintf('bunch charge not given, set to 0 pC, click any key to continue... \n'); q_charge=0.0; %pause;
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% load R56, lattice start and end positions
% format (s,R11,R12,...R16,R21,...,R26,...,R51,R52,R53,R54,R55,R56,...,R65,R66)
delimiterIn=' '; headerlinesIn=0;
transport=importdata(filename03,delimiterIn,headerlinesIn);
s_ele=transport(:,1);                            % Frenet-Serret s in m
R55_ele= transport(:,30);                        % R55, unitless
R56_ele=-transport(:,31);                        % R56, in m
start_pos=s_ele(1);
end_pos=s_ele(end)-0.001;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% evaluate compression factor
%C_ele=1./(1-chirp*R56_ele);
C_ele=1./(R55_ele-chirp*R56_ele);
C_factor=C_ele(end);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% evaluate initial beam current in Amp
if ((isempty(q_charge)==0) && (isempty(sigma_s)==0))
    I_b_bc=q_charge/fwhm_t;
    %I_b=I_b_bc*C_ele(end);
    I_b=I_b_bc;
    fprintf('***calculated initial I_b = %f A, set manually if not desired***\n',I_b);

    disp('set up I_b manually if the above parameter is not desired...');
    
    if (skip_ask==0)
    prompt='press Enter/Return after giving I_b in Amp: (0 if the calculated number is ok...)';
    num=input(prompt);
    if (num ~= 0)
        I_b=num;
    end
    end
else
    fprintf('either bunch charge or sigma_s not given... \n'); I_b=0.0; %pause;
end
%-------------------------------------------------------------------------%

disp('***********************************************************************');
disp('The following summarizes input parameters for running volterra_mat...\n');
fprintf('beam energy in GeV = %f \n',energy);
fprintf('beam current after compression = %f A \n',I_b);
fprintf('final compression factor = %f \n',C_factor);
fprintf('horizontal normalized emittance = %f um \n',emit_norm_x);
fprintf('vertical normalized emittance = %f um \n',emit_norm_y);
fprintf('horizontal initial beta function = %f m \n',betax0);
fprintf('vertical initial beta function = %f m \n',betay0);
fprintf('horizontal initial alpha function = %f \n',alphax0);
fprintf('vertical initial alpha function = %f \n',alphay0);
fprintf('relative energy spread = %f \n',sigma_delta);
fprintf('chirp parameter = %f m^-1 \n',chirp);
fprintf('start position = %f m\n',start_pos);
fprintf('end position = %f m\n',end_pos);
disp('***********************************************************************');

disp('The following simulation parameters should be given. Press Enter to skip to GUI...\n');
fprintf('suggest at least one mesh per 5 cm within a dipole for calculation (suggest to perform convergence test) \n');
fprintf('if CSR drift impedance is included, the minimum mesh size is further limited by formation length \n');
fprintf('total length of the beamline = %f m, shortest dipole length = %f cm...\n',end_pos-start_pos,shortest);
if (shortest > 1.0)
    mesh_num_tmp=3*floor(abs(start_pos-end_pos)*100/shortest);
    fprintf('suggested mesh_num = %d... \n',mesh_num_tmp);
else
    mesh_num_tmp=floor(abs(start_pos-end_pos)*100);
    fprintf('suggested mesh_num = %d... \n',mesh_num_tmp);
end

if (skip_ask==0)
prompt='want to use the suggested mesh_num? (Y:1, N:0)';
tmp=input(prompt);
if (tmp==1)
    mesh_num=mesh_num_tmp;
else
    prompt='give a mesh_num = ';
    mesh_num=input(prompt);
end

fprintf('suggested interval of lambda scan for vanishing emittance: \n');
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Note: this formula only applies for the case of vanishing emittance and 
% for isochrounous lattice. For finite emittance or with bunch compression,
% adding a red-shift for lambda is suggested.
%-------------------------------------------------------------------------%
if (abs(R56_ele(end))<1e-4)
    lambda_opt=3*mean(abs(R56_ele))*sigma_delta*1e6;
else
    lambda_opt=3*abs(R56_ele(end))*sigma_delta*1e6;
end
lambda_start=0.2*lambda_opt;   % unit: um
lambda_end=50*lambda_opt;      % unit: um

fprintf('suggested lambda scan interval: from %f to %f um \n',lambda_start,lambda_end);
%-------------------------------------------------------------------------%
if (skip_ask==0)
fprintf('first interval for lambda scan: \n');
prompt='lambda_start01 (unit: um)=';
lambda_start01=input(prompt);
prompt='lambda_end01 (unit: um)=';
lambda_end01=input(prompt);
prompt='scan_num01= ';
scan_num01=input(prompt);

prompt='second interval for lambda scan: (optional, 0 to suppress)';
opt=input(prompt);
if (opt ~= 0)
    prompt='lambda_start02 (unit: um)=';
    lambda_start02=input(prompt);
    prompt='lambda_end02 (unit: um)=';
    lambda_end02=input(prompt);
    prompt='scan_num02= ';
    scan_num02=input(prompt);
else
    lambda_start02=1;
    lambda_end02=1;
    scan_num02=0;
end

prompt='plot lattice transport function, e.g. R56(s)? (Y:1, N:0)';
iplot_lattice=input(prompt);
prompt='plot beam current evolution I_b(s)? (Y:1, N:0)';
iplot_I_b=input(prompt);
prompt='plot lattice quilt pattern? (Y:1, N:0)';
iquilt_plot=input(prompt);
prompt='plot gain function, i.e. G(s) for a specific lambda? (Y:1, N:0)';
iplot_gain_func=input(prompt);
prompt='plot gain spectrum, i.e. Gf(lambda) at the end of lattice? (Y:1, N:0)';
iplot_gain_spec=input(prompt);
prompt='plot gain map, i.e. G(s,lambda)? (Y:1, N:0)';
surfplot_gain_spec_func=input(prompt);
prompt='plot energy modulation function or spectrum? (Y:1, N:0)';
iplot_energy_mod=input(prompt);
prompt='save as workplace.mat? (Y:1, N:0)';
isave=input(prompt);
prompt='calculate staged (or iterative) solution? (Y:1, N:0)';
iiterative=input(prompt);
if (iiterative==1)
    prompt='calculate stage gain coefficient d_m? (Y:1, N:0)';
    idm_analysis=input(prompt);
else
    idm_analysis=0;
end
prompt='only calculate gain spectrum? (this can speed up calculation) (Y:1, N:0)';
ilast=input(prompt);
prompt='include steady-state CSR in bends? (Y:1, N:0)';
iCSR_ss=input(prompt);
if (iCSR_ss==1)
    prompt='specify ultrarelativistic or non-ultrarelativistic model? (UR:1, NUR:2)';
    ssCSR_model=input(prompt);
    prompt='want to include possible CSR shielding effect? (Y:1, N:0)';
    issCSRpp=input(prompt);
    if (issCSRpp==1)
        prompt='specify the full pipe height in cm...';
        full_pipe_height=input(prompt);
    else
        full_pipe_height=1e50;
    end
else
    ssCSR_model=-1;
    issCSRpp=-1;
end
prompt='include transient CSR in bends? (Y:1, N:0)';
iCSR_tr=input(prompt);
prompt='include CSR in drifts? (Y:1, N:0)';
iCSR_drift=input(prompt);
prompt='include LSC in drifts? (Y:1, N:0)';
iLSC=input(prompt);
if (iLSC==1)
    prompt='specify a model? (1:on-axis, 2:ave, 3:axisymmetric Gaussian)';
    LSC_model=input(prompt);
else
    LSC_model=-1;
end
prompt='include linac geometric impedance? (Y:1, N:0)';
ilinac=input(prompt);
prompt='use coasting or Gaussian model? (1-coasting, 2-Gaussian)';
iGaussian=input(prompt);
fprintf('Gaussian beam model under construction...\n');

prompt='calculate energy modulation function? (Y:1, N:0)';
iEnergy_gain_calc=input(prompt);
prompt='calculate energy modulation spectrum? (Y:1, N:0)';
iEnergy_mod_calc=input(prompt);

prompt='any RF element in the lattice? (Y:1, N:0)';
iRF_ele=input(prompt);
if (iRF_ele==1) 
    trans_mat_trim; 
else
    find_TWLA=0;
end
else
    % default setting
    lambda_start01=1;
    lambda_end01=100;
    scan_num01=10;
    lambda_start02=0;
    lambda_end02=0;
    scan_num02=0;
    mesh_num=600;
    iplot_lattice=0;
    iplot_I_b=0;
    iquilt_plot=0;
    iplot_gain_func=1;
    iplot_gain_spec=1;
    surfplot_gain_spec_func=0;
    iplot_energy_mod=0;
    iCSR_ss=1;
    ssCSR_model=1;
    issCSRpp=0;
    full_pipe_height=1e50;
    round_pipe_radius=1e50;
    iCSR_tr=0;
    iCSR_drift=0;
    iLSC=0;
    LSC_model=1;
    ilinac=0;
    iGaussian=1;
    iEnergy_gain_calc=0;
    iRF_ele=0;
    Derbenev=0;
    iTransverse_gain_calc=0;
    first_harmonic_notification=0;
end
disp('***********************************************************************');

%if (run_mode==1)
    %openfig('/Users/jcytsai/Desktop/matlab-mac/volterra_mat_v1.6/GUI_volterra.fig');
    GUI_volterra;
    
%else
%    volterra_mat_for_direct_run;
%end
