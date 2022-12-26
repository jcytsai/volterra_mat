%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This program calculates microbunching gain based on linearized Vlasov   %
% equation (or, Volterra integral equation in mathematical form), with    %
% collective effects such as coherent synchrotron radiation (CSR),        %
% longitudinal space charge(LSC) and linac geometric impedances considered%
% The impedance models are in analytical form. CSR is assumed 1-D, with   %
% ultra-relativistic (usually used for high-energy case) and relativistic %
% (can be used for low-energy case) versions. There are two options for   %
% CSR models; the free-space model and parallel-plate shielding model. The%
% transient CSR models, including entrance transient and exit transient   %
% (or, CSR drift), are based on Saldin's treatment. The LSC models have   %
% various versions, for different transverse beam distributions and       %
% various shapes (cross sections) of surrounding pipes. See Y.Li's thesis %
% for more details.                                                       %
%                                                                         %
% As of v3.3 and later, we have included both density and energy          %
% induced microbunching. As of v4.0 and later, we have extended theoretical
% treatment to incorporate transverse-longitudinal microbunching. For more%
% detail, they can be found in IPAC16 and NAPAC16 paper.                  %
%                                                                         %
% The required input information includes physical parameters and         %
% numerical parameters. The former are read from ELEGANT, so users need to%
% prepare a set of ELEGANT input files (*.ele and *.lte) for a certain    %
% lattice, and specify a set of ELEGANT output files (*.twi, *.mag, *.sig,%
% *.cen, *.flr, *.mat, *.param, etc) The latter must be provided by users,%
% through the interactive message and GUI.                                %
%                                                                         %
% This program calculates the gain function G(s), if a specific modulation%
% wavelength is given, and/or gain spectrum, if a series of modulation    %
% wavelengths are scanned.                                                %
%                                                                         %
% Since v3.3, we have extended the concepts of microbunching gain function%
% to include:                                                             %
% (i)   density-to-density gain                                           %
% (ii)  density-to-energy gain                                            %
% (iii) energy-to-density gain                                            %
% (iv)  energy-to-energy gain                                             %
% to complete the general microbunching gain analysis.                    %
% This way the definitions of gain functions (also gain spectra) can have %
% four kinds of combinations as well.                                     %
%                                                                         %
% For more details about the original theoretical formulation, please     %
% refer to the paper by Heifets, Stupakov, and Krinsky, PRST-AB 5, 064401 %
% (2002) and the erratum paper. For details about the extended analysis,  %
% pleae refer to our IPAC16 and NAPAC16 work.                             %
%                                                                         %
% Benchmarking of this program against ELEGANT was documented in          %
% JLAB-TN-14-016. There is a step-by-step guide to using this program,    %
% which was documented in JLAB-TN-15-019. Benchmarking against energy     %
% modulation induced microbunching can be found in JLAB-TN-16-022.        %
%                                                                         %
% Details about impedance models used in the code can be found in:        %
% C. -Y. Tsai et al., IPAC'15 (MODMA028)                                  %
%                                                                         %
% Note:                                                                   %
% 1. since v3.3, we have removed iterative (or staged) gain calculation   %
%    from the program.                                                    %
% 2. since v4.0, we have migrated from HSK's formulation to HK's          %
%    formulation.                                                         %
% 3. since v4.0, we have included the microbunching in (z,delta), (x,z),  %
%    (x',z), (y,z), and (y',z).                                           %
%                                                                         %
% Program written by Cheng-Ying Tsai, jcytsai@vt.edu                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
%------------------ Global variable declaration --------------------------%
global egamma s_ele R11_ele R12_ele R13_ele R14_ele R21_ele R22_ele R23_ele R24_ele R31_ele R32_ele R33_ele R34_ele R41_ele R42_ele R43_ele R44_ele R16_ele R26_ele R36_ele R46_ele R51_ele R52_ele R53_ele R54_ele R55_ele R56_ele C_ele;
global rhox rhoy k_wave emit_norm_x emit_norm_y emitx emity alphax0 alphay0 betax0 betay0 gammax0 gammay0 sigma_delta n_1k0 e_1k0 ax_1k0 axp_1k0 ay_1k0 ayp_1k0 re nb I_b chirp start_pos end_pos I_A const_A mesh_num;
%global dipole_s C_factor lambda_start01 lambda_end01 scan_num01 mesh_num;
global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac issCSRpp;
global sig_x_ele sig_y_ele sig_xp_ele sig_yp_ele enx_s_ele eny_s_ele egamma_vec find_TWLA full_pipe_height round_pipe_radius;
global LSC_model ssCSR_model first;
global dipole_control_seq;
global s_ele_Twiss alphax_s_ele alphay_s_ele betax_s_ele betay_s_ele gammax_s_ele gammay_s_ele Nb sigma_z0;
global iIBS s_IBS tau_inv_IBS tau_inv_IBS_x tau_inv_IBS_y emit_gx_IBS emit_gy_IBS sigma_delta_IBS D_IBS_z D_zz;
global DO_X DO_Y DO_Z iDz iDzz itransLD plot_Zk_vs_s iSES CSRdrifmodel Derbenev checkbox31 checkbox32 gen_CSR_LSC_param;

format long
const_A=3^(-1/3)*gamma(2/3)*(1i*sqrt(3)-1);
c_speed=2.99792458e10;
charge=1.6e-19;
first=1; % signal indicating shielding CSR model is used (see kernel_mod.m)
%emit_calc=1; 
%sigma_z0=1e-2;

%iiterative=0;
order_N=3;
%iIBS=0;                      
%DO_X=1;DO_Y=1;DO_Z=1;
%iSES=1; iSES_ana=0;
enhance_factor=1;            % enhance factor for IBS
%iDz=1; iDzz=1;               % 1: on; 0: off
%------------------ Beamline concatenation setting -----------------------%
concatenate=0;
ans_tmp=0;                   %1: write output for next input
file_index_read=2;           %file_index means location index, integer
file_index_write=2;          %remember to remove existing files in advance
                               
if (concatenate==0)          % Note on 20200615: if iIBS=1, n_1k0=1 and the others are 0.
    n_1k0=1.0;
    e_1k0=0.0;
    ax_1k0=0.0;
    axp_1k0=0.0;
    ay_1k0=0.0;
    ayp_1k0=0.0;
end
%------------------ Read input values from GUI ---------------------------%
HH = findobj(gcf,'Tag','energy');          energy = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','I_b');             I_b = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','C_factor');        C_factor = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','emit_norm_x');     emit_norm_x = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','emit_norm_y');     emit_norm_y = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','betax0');          betax0 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','betay0');          betay0 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','alphax0');         alphax0 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','alphay0');         alphay0 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','sigma_delta');     sigma_delta = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','chirp');           chirp = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','start_pos');       start_pos = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','end_pos');         end_pos = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','lambda_start01');  lambda_start01 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','lambda_end01');    lambda_end01 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','scan_num01');      scan_num01 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','lambda_start02');  lambda_start02 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','lambda_end02');    lambda_end02 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','scan_num02');      scan_num02 = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','mesh_num');        mesh_num = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iplot_lattice');   iplot_lattice = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iplot_gain_func'); iplot_gain_func = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iplot_gain_spec'); iplot_gain_spec = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iplot_energy_mod');iplot_energy_mod = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','isave');           isave = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iiterative');      iiterative = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','ilast');           ilast = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iCSR_ss');         iCSR_ss = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iCSR_tr');         iCSR_tr = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iCSR_drift');      iCSR_drift = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iLSC');            iLSC = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','ilinac');          ilinac = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iGaussian');       iGaussian = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iEnergy_mod_calc');iEnergy_mod_calc = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iEnergy_gain_calc');iEnergy_gain_calc = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','idm_analysis');    idm_analysis = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','surfplot_gain_spec_func');surfplot_gain_spec_func = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','issCSRpp');        issCSRpp = str2num(get(HH,'String'));

HH = findobj(gcf,'Tag','iplot_lattice');   iplot_lattice = HH.Value;
HH = findobj(gcf,'Tag','iplot_gain_func'); iplot_gain_func = HH.Value;
HH = findobj(gcf,'Tag','iplot_gain_spec'); iplot_gain_spec = HH.Value;
HH = findobj(gcf,'Tag','iplot_energy_mod');iplot_energy_mod = HH.Value;
HH = findobj(gcf,'Tag','iiterative');      iiterative = HH.Value;
HH = findobj(gcf,'Tag','ilast');           ilast = HH.Value;
HH = findobj(gcf,'Tag','iCSR_ss');         iCSR_ss = HH.Value;
HH = findobj(gcf,'Tag','iCSR_tr');         iCSR_tr = HH.Value;
HH = findobj(gcf,'Tag','iCSR_drift');      iCSR_drift = HH.Value;
HH = findobj(gcf,'Tag','iLSC');            iLSC = HH.Value;
%HH = findobj(gcf,'Tag','ilinac');          ilinac = HH.Value;
HH = findobj(gcf,'Tag','iEnergy_mod_calc');iEnergy_mod_calc = HH.Value;
HH = findobj(gcf,'Tag','iEnergy_gain_calc');iEnergy_gain_calc = HH.Value;
HH = findobj(gcf,'Tag','idm_analysis');    idm_analysis = HH.Value;
HH = findobj(gcf,'Tag','surfplot_gain_spec_func');surfplot_gain_spec_func = HH.Value;
HH = findobj(gcf,'Tag','issCSRpp');        issCSRpp = HH.Value;

HH = findobj(gcf,'Tag','popupmenu6');      ssCSR_model = HH.Value;
HH = findobj(gcf,'Tag','popupmenu7');      CSRdrifmodel = HH.Value;
HH = findobj(gcf,'Tag','popupmenu8');      LSC_model = HH.Value;
HH = findobj(gcf,'Tag','gen_CSR_LSC_param');gen_CSR_LSC_param = HH.Value;
%HH = findobj(gcf,'Tag','ssCSR_model');     ssCSR_model = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','CSRdrifmodel');    CSRdrifmodel = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','LSC_model');       LSC_model = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iquilt_plot');     iquilt_plot = HH.Value;
HH = findobj(gcf,'Tag','iplot_I_b');       iplot_I_b = HH.Value;

HH = findobj(gcf,'Tag','iIBS');            iIBS = HH.Value;
HH = findobj(gcf,'Tag','itransLD');        itransLD = HH.Value;
HH = findobj(gcf,'Tag','plot_Zk_vs_s');    plot_Zk_vs_s = HH.Value;
HH = findobj(gcf,'Tag','iDz');             iDz = HH.Value;
HH = findobj(gcf,'Tag','iDzz');            iDzz = HH.Value;
HH = findobj(gcf,'Tag','iSES');            iSES = HH.Value;
HH = findobj(gcf,'Tag','DO_X');            DO_X = HH.Value;
HH = findobj(gcf,'Tag','DO_Y');            DO_Y = HH.Value;
HH = findobj(gcf,'Tag','DO_Z');            DO_Z = HH.Value;
HH = findobj(gcf,'Tag','Derbenev');        Derbenev = HH.Value;
HH = findobj(gcf,'Tag','checkbox31');      checkbox31 = HH.Value;
HH = findobj(gcf,'Tag','checkbox32');      checkbox32 = HH.Value;
%-------------------------------------------------------------------------%

save_SDDS=checkbox32;

iSES_ana=iSES;
%if (ilinac==1); iRF_ele=1; end
ilinac=0; iRF_ele=0;

%HH = findobj(gcf,'Tag','ssCSR_model');     ssCSR_model = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','LSC_model');       LSC_model = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','round_pipe_radius');round_pipe_radius = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iquilt_plot');     iquilt_plot = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iplot_I_b');       iplot_I_b = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iRF_ele');         iRF_ele = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','full_pipe_height');full_pipe_height = str2num(get(HH,'String'));

iTransverse_gain_calc=0;
first_harmonic_notification=0;

%HH = findobj(gcf,'Tag','Derbenev');        Derbenev = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iTransverse_gain_calc');iTransverse_gain_calc = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','first_harmonic_notification');first_harmonic_notification = str2num(get(HH,'String'));

if (iIBS==1 && DO_X==0 && DO_Y==0 && DO_Z==0)
    fprintf('something wrong in IBS setting...\n'); pause;
end

if (iIBS==1)
    iEnergy_gain_calc=1;
    iTransverse_gain_calc=0;             % set zero due to under construction
    iiterative=0;                        % required
    D_friction=1;
    D_diffusion=1;
end

if (iTransverse_gain_calc==1 && iEnergy_gain_calc==0)
    fprintf('to calculate transverse gain, energy gain must be calculated. iEnergy_gain_calc automatically set to 1.\n');
    iEnergy_gain_calc=1;
end
%{
%--------------------- Initialize usage of mtimesx -----------------------%
blas_lib='/Applications/MATLAB_R2015b.app/bin/maci64/libmwblas.dylib';
mex('-DDEFINEUNIX','-largeArrayDims','mtimesx.c',blas_lib);
clc;
%}

%--------------------- Initialize beam parameters ------------------------%
egamma=energy*1000/0.511;                           % Lorentz factor
I_A=17045;                                          % Alfven current in Amp
re=2.81794*10^(-13);                                % classical radius in cm
nb=I_b/(re*I_A);                                    % line density (cm^-1)
Nb=2*sqrt(2*log(2))*sigma_z0*I_b/(c_speed*charge);  % total number of physical particles
emit_norm_x=emit_norm_x*10^(-4);                    % normalized horizontal emittance in cm
emitx=emit_norm_x/egamma;                           % geometric  horizontal emittance in cm
emit_norm_y=emit_norm_y*10^(-4);                    % normalized vertical   emittance in cm
emity=emit_norm_y/egamma;                           % geometric  vertical   emittance in cm
betax0=betax0*10^2;                                 % initial betax in cm
betay0=betay0*10^2;                                 % initial betay in cm
gammax0=(1+alphax0^2)/betax0;                       % initial gammax in cm^(-1)
gammay0=(1+alphay0^2)/betay0;                       % initial gammay in cm^(-1)
chirp=chirp/100;                                    % chirp parameter in cm^(-1), negative for compression

%--------------------- Initialize lattice parameters ---------------------%
start_pos=start_pos*10^2;
end_pos=end_pos*10^2;
rhox=rhox*10^2;                                     % horizontal bending radius in cm
rhoy=rhoy*10^2;                                     % vertical bending radius in cm

call_dipole;
load('dipole_control_seq.txt');

[r1,~]=size(dipole_control_seq);
[r2,~]=size(dipole_s);

if (r1~=(r2/2))
    fprintf('something inconsistent with setting of dipole_control_seq...\n');
    fprintf('ignore the message if care of dipole action control is not taken...\n');
end

load('lattice_transport_functions.o');
transport=lattice_transport_functions;
load('pcentral_function.o');
egamma_vec=interp1(pcentral_function(:,1),pcentral_function(:,2),transport(:,1));

% format (s,R11,R12,...,R51,R52,R53,R54,R55,R56,...,R65,R66) with 1+36 columes of data
% use coordinates (x,px,y,py,z,delta_gamma)
% delimiterIn='\t';
delimiterIn=' '; 
headerlinesIn=0;
%transport=importdata(filename,delimiterIn,headerlinesIn);
    
s_ele=100*transport(:,1);                           % Frenet-Serret s (in cm)
R11_ele=transport(:,2);
R12_ele=100*transport(:,3);
R13_ele=transport(:,4);
R14_ele=100*transport(:,5);
R16_ele=100*transport(:,7);                         % dispersion Dx(s)
R21_ele=transport(:,8)/100;
R22_ele=transport(:,9);
R23_ele=transport(:,10)/100;
R24_ele=transport(:,11);
R26_ele=transport(:,13);
R31_ele=transport(:,14);
R32_ele=100*transport(:,15);
R33_ele=transport(:,16);
R34_ele=100*transport(:,17);
R36_ele=100*transport(:,19);                        % dispersion Dy(s)
R41_ele=transport(:,20)/100;
R42_ele=transport(:,21);
R43_ele=transport(:,22)/100;
R44_ele=transport(:,23);
R46_ele=transport(:,25);
R51_ele=transport(:,26);                            % R51 unitless
R52_ele=100*transport(:,27);                        % R52 in cm
R53_ele=transport(:,28);                            % R53 unitless
R54_ele=100*transport(:,29);                        % R54 in cm
R55_ele=transport(:,30);                            % R55 unitless
R56_ele=-100*transport(:,31);                       % R56 in cm
%{
%********* conversion of R51: from ELEGANT to HSK format *****************%
R51_ele=-R51_ele+R52_ele*(alphax0/betax0);
R52_ele=-R52_ele;
if (betay0 > 0)
    R53_ele=-R53_ele+R54_ele*(alphay0/betay0);
    R54_ele=-R54_ele;
end
%*************************************************************************%
%}
%********* conversion of R5x: from ELEGANT to HK format ******************%
R51_ele=-R51_ele;
R52_ele=-R52_ele;
R53_ele=-R53_ele;
R54_ele=-R54_ele;
%*************************************************************************%
%
%******** additional contribution of R56 due to drfit section ************%
%if ((iCSR_drift==1)||(iLSC==1)||(ilinac==1))
	R56_drift=(s_ele-start_pos)./egamma_vec.^2;
	R56_ele=R56_ele+R56_drift;
%end
%*************************************************************************%
%}

%*************************************************************************%
if (iRF_ele==0); find_TWLA=0; end

C_ele=1./(R55_ele-chirp*R56_ele);
C_ele=abs(C_ele);                                   % avoid roll-over or parasitic compression

%******* Read transverse rms beam size function for LSC calculation ******%
filename='beam_sigma_mat_functions.o';
% foramt (s,sigma_x,sigma_y,sigma_s,sigma_delta,emit_norm_x,emit_norm_y)
beam_sigma_mat_tmp=importdata(filename,delimiterIn,headerlinesIn);    
beam_sigma_mat=interp1(beam_sigma_mat_tmp(:,1),beam_sigma_mat_tmp,s_ele/100);
sig_x_ele=100*beam_sigma_mat(:,2);                  % sigma_x in cm
sig_y_ele=100*beam_sigma_mat(:,3);                  % sigma_y in cm
sig_s_ele=100*beam_sigma_mat(:,4);                  % sigma_s in cm
%sdelta_s_ele= beam_sigma_mat(:,5);
enx_s_ele=100*beam_sigma_mat(:,6);                  % enx_s in cm
eny_s_ele=100*beam_sigma_mat(:,7);                  % eny_s in cm

sig_xp_ele=beam_sigma_mat(:,8);                  % sigma_xp in rad
sig_yp_ele=beam_sigma_mat(:,9);                  % sigma_yp in rad

C_ele_alt=sig_s_ele(1)./sqrt(R56_ele.^2*sigma_delta^2+(R55_ele-chirp*R56_ele).^2*(sig_s_ele(1))^2);
C_ele=C_ele_alt;

%C_ele_sigs=sig_s_ele(1)./sig_s_ele;
%C_ele=C_ele_sigs;
%C_ele=sig_s_ele(1)./sig_s_ele;

t=0; while (isnan(enx_s_ele(end-t))); t=t+1; end
if (t>0); enx_s_ele(end-t+1:end,1)=ones(t,1)*enx_s_ele(end-t,1); end
t=0; while (isnan(eny_s_ele(end-t))); t=t+1; end
if (t>0); eny_s_ele(end-t+1:end,1)=ones(t,1)*eny_s_ele(end-t,1); end
t=0; while (isnan(sig_x_ele(end-t))); t=t+1; end
if (t>0); sig_x_ele(end-t+1:end,1)=ones(t,1)*sig_x_ele(end-t,1); end
t=0; while (isnan(sig_y_ele(end-t))); t=t+1; end
if (t>0); sig_y_ele(end-t+1:end,1)=ones(t,1)*sig_y_ele(end-t,1); end
t=0; while (isnan(sig_xp_ele(end-t))); t=t+1; end
if (t>0); sig_xp_ele(end-t+1:end,1)=ones(t,1)*sig_xp_ele(end-t,1); end
t=0; while (isnan(sig_yp_ele(end-t))); t=t+1; end
if (t>0); sig_yp_ele(end-t+1:end,1)=ones(t,1)*sig_yp_ele(end-t,1); end
t=0; while (isnan(sig_s_ele(end-t))); t=t+1; end
if (t>0); sig_s_ele(end-t+1:end,1)=ones(t,1)*sig_s_ele(end-t,1); end
%*************************************************************************%

%************************** Plot lattice transport functions *************%
if (iplot_lattice==1)
	figure(101); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R16_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{16} (cm)');       grid off; hold on; axis('tight');
	figure(102); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R36_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{36} (cm)');       grid off; hold on; axis('tight');
	figure(103); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R51_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{51}');            grid off; hold on; axis('tight');
	figure(104); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R52_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{52} (cm)');       grid off; hold on; axis('tight');
	figure(105); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R53_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{53}');            grid off; hold on; axis('tight');
	figure(106); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R54_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{54} (cm)');       grid off; hold on; axis('tight');
    figure(107); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R56_ele,'b-','linewidth',5);   xlabel('s (m)'); ylabel('R_{56} (cm)');       grid off; hold on; axis('tight');
    figure(108); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele,'b-','linewidth',5);     xlabel('s (m)'); ylabel('Compression factor');grid off; hold on; axis('tight');
end
%*************************************************************************%

%************************** Plot beam current evolution ******************%
if (iplot_I_b==1)
    figure(109); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele*I_b,'b-','linewidth',5); xlabel('s (m)'); ylabel('Peak current (A)'); grid off; hold on; axis('tight'); 
end
%*************************************************************************%

%
%******* Read Twiss functions for emittance calculation ******************%
% UPDATE on 20200224: below used for IBS calculation, see g0k_mat.m
%if (emit_calc==1)
    filename02='Twiss_lattice.o';
    % foramt (s,alphax,alphay,betax,betay,psix,psiy)
    Twiss_mat_ele=importdata(filename02,delimiterIn,headerlinesIn);
    
    [r,~]=size(Twiss_mat_ele);
    p=2; Twiss_mat(1,:)=Twiss_mat_ele(1,:);
    for m=2:1:r
        if (Twiss_mat_ele(m,1)>Twiss_mat_ele(m-1,1))
            Twiss_mat(p,:)=Twiss_mat_ele(m,:);
            p=p+1;
        end
    end
    
    s_ele_Twiss=Twiss_mat(:,1)*1e2;                 % unit: cm
    alphax_s_ele=Twiss_mat(:,2);
    alphay_s_ele=Twiss_mat(:,3);
    betax_s_ele=Twiss_mat(:,4)*1e2;                 % unit: cm
    betay_s_ele=Twiss_mat(:,5)*1e2;                 % unit: cm
    gammax_s_ele=(1+alphax_s_ele.^2)./betax_s_ele;  % unit: cm^-1
    gammay_s_ele=(1+alphax_s_ele.^2)./betax_s_ele;  % unit: cm^-1    
    %s=linspace(start_pos,end_pos,mesh_num);
    %xemit_calc(s);
    %yemit_calc(s); 
%end
%}
%************************** Plot quilt pattern R56(s'->s) ****************%
if (iquilt_plot==1) 
    quilt_plot; 
end
%*************************************************************************%

%--------------------- Initialize numerical parameters -------------------%
lambda_start01=lambda_start01*10^(-4);              % starting modulation wavelength in cm
lambda_end01  =lambda_end01*10^(-4);                % ending modulation wavelength in cm                                         
scan_mesh01   =(lambda_end01-lambda_start01)/scan_num01;  
lambda_start02=lambda_start02*10^(-4);              % starting modulation wavelength in cm
lambda_end02  =lambda_end02*10^(-4);                % ending modulation wavelength in cm
scan_mesh02   =(lambda_end02-lambda_start02)/scan_num02;
scan_num      =scan_num01+scan_num02;
mesh_size     =(end_pos-start_pos)/mesh_num;        % integration steps in cm, uniform mesh

%--------------------- Initialize progress bar ---------------------------%
%WB=waitbar(0,'Please wait...');
progressbar; 
% This subroutine is written by Steve Hoelzer, see
% http://www.mathworks.com/matlabcentral/fileexchange/6922-progressbar

%--------------------- IBS calculation -----------------------------------%
    if (iIBS==1)
        s=linspace(start_pos,end_pos,mesh_num);
        s_IBS=s; 
        %IBS_calc;                           % this way may lead to error when refresh memory
        %
        use_Piwinski=0;                      % Piwinski model
        use_CIMP=1;                          % Bjorken-Mtingwa model
        mod_factor=2*sqrt(2*log(2));         % for Gaussian bunch
        %mod_factor=1;                       % for uniform bunch
        update_SES=1;                        % update slice energy spread
        
        scale_factor=1;
        tmp01=interp1(s_ele,C_ele,s);               % C(s)
        ds=s(2)-s(1);
        sigma_delta_IBS=sigma_delta*ones(1,length(s));
        emit_gx_IBS=emitx*ones(1,length(s));
        emit_gy_IBS=emity*ones(1,length(s));
        C_array=interp1(s_ele,C_ele,s_IBS);
        dC_array=gradient(C_array)./C_array;
        
        if (use_CIMP==1)
            const_factor_01=re^2/(c_speed*charge);
        
            for m=1:1:length(s)
                if (m==1)
                    sigma_delta_IBS(1,1)=sigma_delta;
                    emit_gx_IBS(1,1)=emitx;
                    emit_gy_IBS(1,1)=emity;
                else
                    if (DO_Z==1); sigma_delta_IBS(1,m)=sigma_delta_IBS(1,m-1)*exp(real(tau_inv_IBS(m-1))*ds); end
                    if (DO_X==1); emit_gx_IBS(1,m)=emit_gx_IBS(1,m-1)*exp(2*real(tau_inv_IBS_x(m-1))*ds); end
                    if (DO_Y==1); emit_gy_IBS(1,m)=emit_gy_IBS(1,m-1)*exp(2*real(tau_inv_IBS_y(m-1))*ds); end
                end
                
                s_loc=s(m);
                egamma_tmp=interp1(s_ele,egamma_vec,s_loc);
                betax=interp1(s_ele_Twiss,betax_s_ele,s_loc);
                betay=interp1(s_ele_Twiss,betay_s_ele,s_loc);
                alphax=interp1(s_ele_Twiss,alphax_s_ele,s_loc);
                alphay=interp1(s_ele_Twiss,alphay_s_ele,s_loc);
                R16=interp1(s_ele,R16_ele,s_loc);
                R26=interp1(s_ele,R26_ele,s_loc);
                R36=interp1(s_ele,R36_ele,s_loc);
                R46=interp1(s_ele,R46_ele,s_loc);
                sigx=interp1(s_ele,sig_x_ele,s_loc,'linear','extrap');
                sigy=interp1(s_ele,sig_y_ele,s_loc,'linear','extrap');
                emit_nx=interp1(s_ele,enx_s_ele,s_loc,'linear','extrap'); emit_gx=emit_nx/egamma_tmp;
                emit_ny=interp1(s_ele,eny_s_ele,s_loc,'linear','extrap'); emit_gy=emit_ny/egamma_tmp;
                
                tmp12=1/mod_factor/(64*pi^2);
                
                %A_coeff(m)=tmp12*const_factor_01*(I_b*tmp01(m))/(egamma_tmp^2*emit_nx*emit_ny*sigma_delta_IBS(1,m));
                A_coeff(m)=tmp12*const_factor_01*(I_b*tmp01(m))/(egamma_tmp^2*emit_gx_IBS(1,m)*egamma_tmp*emit_gy_IBS(1,m)*egamma_tmp*sigma_delta_IBS(1,m));
                
                %log_factor(m)=log(egamma_tmp^2*sigy*emit_gx/(re*betax));
                log_factor(m)=log(egamma_tmp^2*sigy*emit_gx_IBS(1,m)/(re*betax));
                H_x(m)=(R16^2+(betax*R26+alphax*R16)^2)/betax;
                H_y(m)=(R36^2+(betay*R46+alphay*R36)^2)/betay;
                %sigma_H(m)=sqrt(1/(1/sigma_delta_IBS(1,m)^2+H_x(m)/emit_gx+H_y(m)/emit_gy));
                sigma_H(m)=sqrt(1/(1/sigma_delta_IBS(1,m)^2+H_x(m)/emit_gx_IBS(1,m)+H_y(m)/emit_gy_IBS(1,m)));
                %a_param(m)=sigma_H(m)/egamma_tmp*sqrt(betax/emit_gx);
                %b_param(m)=sigma_H(m)/egamma_tmp*sqrt(betay/emit_gy);
                a_param(m)=sigma_H(m)/egamma_tmp*sqrt(betax/emit_gx_IBS(1,m));
                b_param(m)=sigma_H(m)/egamma_tmp*sqrt(betay/emit_gy_IBS(1,m));
                
                %{
                % below are Eqs.(38) to (40) of Kubo,Mtingwa,Wolski, PRST-AB 8, 081001 (2005), valid for flat beam only
                g_func_ba=Piwinski_g_func(b_param(m)/a_param(m));
                g_func_ab=Piwinski_g_func(a_param(m)/b_param(m));
                g_ratio_sum=g_func_ba/a_param(m)+g_func_ab/b_param(m);
                tau_inv_IBS(m)=enhance_factor*2*2*pi^(3/2)*A_coeff(m)*log_factor(m)*(sigma_H(m)^2/sigma_delta_IBS(1,m)^2*g_ratio_sum);
                tau_inv_IBS(m)=tau_inv_IBS(m)+dC_array(m)/ds;  % include bunch compression, a factor of 2 ahead for lack of synchrotron motion
                tau_inv_IBS_x(m)=enhance_factor*2*pi^(3/2)*A_coeff(m)*log_factor(m)*(-a_param(m)*g_func_ba+H_x(m)*sigma_H(m)^2/emit_gx*g_ratio_sum);
                tau_inv_IBS_y(m)=enhance_factor*2*pi^(3/2)*A_coeff(m)*log_factor(m)*(-b_param(m)*g_func_ab+H_y(m)*sigma_H(m)^2/emit_gy*g_ratio_sum);
                %}
                %
                % below are Eqs.(35) to (37) of Kubo,Mtingwa,Wolski, PRST-AB 8, 081001 (2005), valid for flat&round beam
                log_factor_01(m)=log(egamma_tmp^2*min(sigx,sigy)*emit_gx_IBS(1,m)/(re*betax));
                log_factor_02(m)=log(egamma_tmp^2*min(sigx,sigy)*emit_gy_IBS(1,m)/(re*betay));
                g_func_ba=log_factor_01(m)*Piwinski_g_func(b_param(m)/a_param(m));
                g_func_ab=log_factor_02(m)*Piwinski_g_func(a_param(m)/b_param(m));
                g_ratio_sum=g_func_ba/a_param(m)+g_func_ab/b_param(m);
                tau_inv_IBS(m)=enhance_factor*2*2*pi^(3/2)*A_coeff(m)*(sigma_H(m)^2/sigma_delta_IBS(1,m)^2*g_ratio_sum);
                tau_inv_IBS(m)=tau_inv_IBS(m)+dC_array(m)/ds;  % include bunch compression, a factor of 2 ahead for lack of synchrotron motion
                tau_inv_IBS_x(m)=enhance_factor*2*pi^(3/2)*A_coeff(m)*(-a_param(m)*g_func_ba+H_x(m)*sigma_H(m)^2/emit_gx*g_ratio_sum);
                tau_inv_IBS_y(m)=enhance_factor*2*pi^(3/2)*A_coeff(m)*(-b_param(m)*g_func_ab+H_y(m)*sigma_H(m)^2/emit_gy*g_ratio_sum);
                %}
            end
        end
                
        egamma_tmp=interp1(s_ele,egamma_vec,s);
        CLog=log_factor;                                        % vector
        const_factor_01=sqrt(pi)*CLog*re./(4*egamma_tmp.^2);    % vector
        Lambda_factor=re*nb./egamma_tmp.*tmp01;                 % vector
        sigx_IBS=interp1(s_ele,sig_x_ele,s_IBS,'linear','extrap');
        sigy_IBS=interp1(s_ele,sig_y_ele,s_IBS,'linear','extrap');
        sigxp_IBS=interp1(s_ele,sig_xp_ele,s_IBS,'linear','extrap');
        sigyp_IBS=interp1(s_ele,sig_yp_ele,s_IBS,'linear','extrap');
        %D_Stupakov=1*const_factor_01*Lambda_factor./(sqrt(sigxp_IBS.*sigyp_IBS).*sigx_IBS.*sigy_IBS);
        
        %
        % transverse Gaussian
        D_zz=iDzz*enhance_factor*const_factor_01.*Lambda_factor./(sqrt(sigxp_IBS.*sigyp_IBS).*sigx_IBS.*sigy_IBS);
        D_IBS_z=iDz*2*enhance_factor*CLog*re^2*nb.*tmp01./(egamma^4*sigx_IBS.*sigy_IBS.*sigxp_IBS.*sigyp_IBS);
        
    end

%--------------------- Option for concatenation study --------------------%
if (concatenate==1)
    if (iEnergy_gain_calc==1 && iTransverse_gain_calc==0)
        
        [tmp01,tmp02,tmp03,tmp04]=read_tot_col_vec_4x4(2,scan_num,file_index_read);
        n_1k0_array=(tmp01+tmp02);
        e_1k0_array=(tmp03+tmp04);
        
    elseif (iEnergy_gain_calc==1 && iTransverse_gain_calc==1)
        
        [tmp01,tmp02,tmp03,tmp04,tmp05,tmp06,tmp07,tmp08,tmp09,tmp10,tmp11,tmp12,...
         tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,tmp20,tmp21,tmp22,tmp23,tmp24,...
         tmp25,tmp26,tmp27,tmp28,tmp29,tmp30,tmp31,tmp32,tmp33,tmp34,tmp35,tmp36]=read_tot_col_vec_36x36(2,scan_num,file_index_read);
     
        n_1k0_array=(tmp01+tmp02+tmp03+tmp04+tmp05+tmp06);
        e_1k0_array=(tmp07+tmp08+tmp09+tmp10+tmp11+tmp12);
        ax_1k0_array=(tmp13+tmp14+tmp15+tmp16+tmp17+tmp18);
        axp_1k0_array=(tmp19+tmp20+tmp21+tmp22+tmp23+tmp24);
        ay_1k0_array=(tmp25+tmp26+tmp27+tmp28+tmp29+tmp30);
        ayp_1k0_array=(tmp31+tmp32+tmp33+tmp34+tmp35+tmp36);
        
    end
    %lambda_array=linspace(lambda_start01,lambda_end01,scan_num);
    %figure(30); [hAx,hLine1,hLine2]=plotyy(lambda_array*1e4,abs(n_1k0_array),lambda_array*1e4,abs(e_1k0_array)); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('\lambda (\mum)'); ylabel(hAx(1),'n_{1,k}(0)'); ylabel(hAx(2),'e_{1,k}(0)');
    %figure(31); [hAx,hLine1,hLine2]=plotyy(lambda_array*1e4,real(n_1k0_array),lambda_array*1e4,imag(n_1k0_array)); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('\lambda (\mum)'); ylabel(hAx(1),'Re n_{1,k}(0)'); ylabel(hAx(2),'Im n_{1,k}(0)');
    %figure(32); [hAx,hLine1,hLine2]=plotyy(lambda_array*1e4,real(e_1k0_array),lambda_array*1e4,imag(e_1k0_array)); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('\lambda (\mum)'); ylabel(hAx(1),'Re e_{1,k}(0)'); ylabel(hAx(2),'Im e_{1,k}(0)');
end

%--------------------- Start of the main program -------------------------%
for iscan=1:1:scan_num
    
    if (scan_num==1)                                      
        lambda=lambda_start01;
    else
        if (iscan <= scan_num01)
            lambda=lambda_start01+(iscan-1)*scan_mesh01;
        else
            lambda=lambda_start02+(iscan-1-scan_num01)*scan_mesh02;
        end
    end
    
    progressbar(iscan/scan_num);
    
    lambda_array(iscan)=lambda;
    k_wave=2*pi/lambda;                             % modulation wave number
    
    %{
    if (concatenate==0)
        n_1k0=1.0;                                  % a constant
        e_1k0=0.0;                                  % a constant
        if (iTransverse_gain_calc==1)
            ax_1k0=0.0;
            axp_1k0=0.0;
            ay_1k0=0.0;
            ayp_1k0=0.0;
        end
    else
    %}
    if (concatenate==1)
        n_1k0=n_1k0_array(iscan);
        e_1k0=e_1k0_array(iscan);
        if (iTransverse_gain_calc==1)
            ax_1k0=ax_1k0_array(iscan);
            axp_1k0=axp_1k0_array(iscan);
            ay_1k0=ay_1k0_array(iscan);
            ayp_1k0=ayp_1k0_array(iscan);
        end
    end
%----------------------Solve Volterra integral equation-------------------%
    
    s=linspace(start_pos,end_pos,mesh_num);
    %{
    if (iIBS==1)
        s_IBS=s; 
        IBS_calc; 
    end
    %}
    G0_k=g0k_mat(s);                                % density bunching due to density modulation
    if (iEnergy_gain_calc==1)
        G0_kp=g0kp_mat(s);                          % density bunching due to energy  modulation
        E0_k=p0k_mat(s);                            % energy  bunching due to density modulation
        E0_kp=p0kp_mat(s);                          % energy  bunching due to energy  modulation
        
        if (iTransverse_gain_calc==1)
            G0_kxz=g0kxz_mat(s);                    % density bunching due to (x,z)  modulation
            G0_kxpz=g0kxpz_mat(s);                  % density bunching due to (x',z) modulation
            G0_kyz=g0kyz_mat(s);                    % density bunching due to (y,z)  modulation
            G0_kypz=g0kypz_mat(s);                  % density bunching due to (y',z) modulation
            
            E0_kxz=p0kxz_mat(s);                    % energy  bunching due to (x,z)  modulation
            E0_kxpz=p0kxpz_mat(s);                  % energy  bunching due to (x',z) modulation
            E0_kyz=p0kyz_mat(s);                    % energy  bunching due to (y,z)  modulation
            E0_kypz=p0kypz_mat(s);                  % energy  bunching due to (y',z) modulation
            
            Ax0_k=ax0k_mat(s);                      % bunching in (x,z) due to density modulation
            Ax0_kp=ax0kp_mat(s);                    % bunching in (x,z) due to energy  modulation
            Ax0_kxz=ax0kxz_mat(s);                  % bunching in (x,z) due to (x,z)   modulation
            Ax0_kxpz=ax0kxpz_mat(s);                % bunching in (x,z) due to (x',z)  modulation
            Ax0_kyz=ax0kyz_mat(s);                  % bunching in (x,z) due to (y,z)   modulation
            Ax0_kypz=ax0kypz_mat(s);                % bunching in (x,z) due to (y',z)  modulation
            
            Axp0_k=axp0k_mat(s);                    % bunching in (x',z) due to density modulation
            Axp0_kp=axp0kp_mat(s);                  % bunching in (x',z) due to energy  modulation
            Axp0_kxz=axp0kxz_mat(s);                % bunching in (x',z) due to (x,z)   modulation
            Axp0_kxpz=axp0kxpz_mat(s);              % bunching in (x',z) due to (x',z)  modulation
            Axp0_kyz=axp0kyz_mat(s);                % bunching in (x',z) due to (y,z)   modulation
            Axp0_kypz=axp0kypz_mat(s);              % bunching in (x',z) due to (y',z)  modulation
            
            Ay0_k=ay0k_mat(s);                      % bunching in (y,z) due to density modulation
            Ay0_kp=ay0kp_mat(s);                    % bunching in (y,z) due to energy  modulation
            Ay0_kxz=ay0kxz_mat(s);                  % bunching in (y,z) due to (x,z)   modulation
            Ay0_kxpz=ay0kxpz_mat(s);                % bunching in (y,z) due to (x',z)  modulation
            Ay0_kyz=ay0kyz_mat(s);                  % bunching in (y,z) due to (y,z)   modulation
            Ay0_kypz=ay0kypz_mat(s);                % bunching in (y,z) due to (y',z)  modulation
            
            Ayp0_k=ayp0k_mat(s);                    % bunching in (y',z) due to density modulation
            Ayp0_kp=ayp0kp_mat(s);                  % bunching in (y',z) due to energy  modulation
            Ayp0_kxz=ayp0kxz_mat(s);                % bunching in (y',z) due to (x,z)   modulation
            Ayp0_kxpz=ayp0kxpz_mat(s);              % bunching in (y',z) due to (x',z)  modulation
            Ayp0_kyz=ayp0kyz_mat(s);                % bunching in (y',z) due to (y,z)   modulation
            Ayp0_kypz=ayp0kypz_mat(s);              % bunching in (y',z) due to (y',z)  modulation            
        end
    end
    
    q=1;                                            % dipole index
    K_mat=zeros(mesh_num);
    
    if (iEnergy_gain_calc==1 && iTransverse_gain_calc==0)
        M_minus_L_mat=zeros(mesh_num);
        K_mat_IBS_zz_1=zeros(mesh_num);
        K_mat_IBS_zz_2=zeros(mesh_num);
        K_mat_IBS_zz_30=zeros(mesh_num);            % without sigma_delta_IBS attached
        K_mat_IBS_zz_31=zeros(mesh_num);            % with sigma_delta_IBS attached
        K_mat_IBS_zz_4=zeros(mesh_num);
        K_mat_IBS_z_0=zeros(mesh_num);
        K_mat_IBS_z_1=zeros(mesh_num);
        K_mat_IBS_z_0_pp=zeros(mesh_num);
        K_mat_IBS_z_1_pp=zeros(mesh_num);
        K_mat_IBS_z_2_pp=zeros(mesh_num);
    elseif (iEnergy_gain_calc==1 && iTransverse_gain_calc==1)
        A_mat=zeros(mesh_num);
        B_mat=zeros(mesh_num);
        C_mat=zeros(mesh_num);
        D_mat=zeros(mesh_num);
        MBI_kernel_mat=eye(36*mesh_num,36*mesh_num);
    end
    
    if (Derbenev==1)
        for p=1:1:(mesh_num-1)
            Derbenev_ratio(iscan,p)=interp1(s_ele,sig_x_ele,s(p))/(lambda^(2/3)*abs(auxr(s(p)))^(1/3));
        end
    end
    
    if ((iCSR_drift==1) || (iLSC==1) || (ilinac==1))
    for p=1:1:(mesh_num-1)  % p: column index
    	s_q=s(p+1:mesh_num);
        if (p==1)
            K_mat(p+1:mesh_num,p)=0.5*kernel_mod(s_q,s(p));  
            if (iEnergy_gain_calc==1)
                M_minus_L_mat(p+1:mesh_num,p)=0.5*(kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p)));
                %M_minus_L_mat(p+1:mesh_num,p)=0.5*(-kernel_mod_L(s_q,s(p)));
                if (iTransverse_gain_calc==1)
                    A_mat(p+1:mesh_num,p)=0.5*kernel_mod_A(s_q,s(p));
                    B_mat(p+1:mesh_num,p)=0.5*kernel_mod_B(s_q,s(p));
                    C_mat(p+1:mesh_num,p)=0.5*kernel_mod_C(s_q,s(p));
                    D_mat(p+1:mesh_num,p)=0.5*kernel_mod_D(s_q,s(p));
                end
            end
        else
            K_mat(p+1:mesh_num,p)=kernel_mod(s_q,s(p));
            if (iEnergy_gain_calc==1)
                M_minus_L_mat(p+1:mesh_num,p)=kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p));
                %M_minus_L_mat(p+1:mesh_num,p)=-kernel_mod_L(s_q,s(p));
                if (iTransverse_gain_calc==1)
                    A_mat(p+1:mesh_num,p)=kernel_mod_A(s_q,s(p));
                    B_mat(p+1:mesh_num,p)=kernel_mod_B(s_q,s(p));
                    C_mat(p+1:mesh_num,p)=kernel_mod_C(s_q,s(p));
                    D_mat(p+1:mesh_num,p)=kernel_mod_D(s_q,s(p));
                end
            end
        end
    end
    
    else
    
    % below assume CSR occurs only within dipoles so filling K_mat can further speed    
    for p=1:1:(mesh_num-1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (auxr(s(p)) > 10^20)   % If auxr does not locate within dipoles, the CSR kernel vanishes.
                                  % This condition is not required, but it is for speeding up calculation.
            K_mat(p+1:mesh_num,p)=0.0;
            if (iEnergy_gain_calc==1)
                M_minus_L_mat(p+1:mesh_num,p)=0.0;
                if (iTransverse_gain_calc==1)
                    A_mat(p+1:mesh_num,p)=0.0;
                    B_mat(p+1:mesh_num,p)=0.0;
                    C_mat(p+1:mesh_num,p)=0.0;
                    D_mat(p+1:mesh_num,p)=0.0;
                end
            end
            q=1;
        else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            p_dipole(q)=p; q=q+1;
            s_q=s(p+1:mesh_num);            
            if (p==1)
                K_mat(p+1:mesh_num,p)=0.5*kernel_mod(s_q,s(p));
                if (iEnergy_gain_calc==1)
                    M_minus_L_mat(p+1:mesh_num,p)=0.5*(kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p)));
                    %M_minus_L_mat(p+1:mesh_num,p)=0.5*(-kernel_mod_L(s_q,s(p)));
                    if (iTransverse_gain_calc==1)
                        A_mat(p+1:mesh_num,p)=0.5*kernel_mod_A(s_q,s(p));
                        B_mat(p+1:mesh_num,p)=0.5*kernel_mod_B(s_q,s(p));
                        C_mat(p+1:mesh_num,p)=0.5*kernel_mod_C(s_q,s(p));
                        D_mat(p+1:mesh_num,p)=0.5*kernel_mod_D(s_q,s(p));
                    end
                end
            else
                K_mat(p+1:mesh_num,p)=kernel_mod(s_q,s(p));
                if (iEnergy_gain_calc==1)
                    M_minus_L_mat(p+1:mesh_num,p)=kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p));
                    %M_minus_L_mat(p+1:mesh_num,p)=-kernel_mod_L(s_q,s(p));
                    if (iTransverse_gain_calc==1)
                        A_mat(p+1:mesh_num,p)=kernel_mod_A(s_q,s(p));
                        B_mat(p+1:mesh_num,p)=kernel_mod_B(s_q,s(p));
                        C_mat(p+1:mesh_num,p)=kernel_mod_C(s_q,s(p));
                        D_mat(p+1:mesh_num,p)=kernel_mod_D(s_q,s(p));
                    end
                end
            end
        end
        
    end
    end
    %
    if (iIBS==1)
    for p=1:1:(mesh_num-1)  % p: column index
    	s_q=s(p+1:mesh_num);
        if (p==1)
            if (D_diffusion==1)
            K_mat_IBS_zz_1(p+1:mesh_num,p)  =0.5*kernel_mod_IBS_zz(s_q,s(p),1);
            K_mat_IBS_zz_2(p+1:mesh_num,p)  =0.5*kernel_mod_IBS_zz(s_q,s(p),2);
            K_mat_IBS_zz_30(p+1:mesh_num,p) =0.5*kernel_mod_IBS_zz(s_q,s(p),30);
            K_mat_IBS_zz_31(p+1:mesh_num,p) =0.5*kernel_mod_IBS_zz(s_q,s(p),31);
            K_mat_IBS_zz_4(p+1:mesh_num,p)  =0.5*kernel_mod_IBS_zz(s_q,s(p),4);
            end
            if (D_friction==1)
            K_mat_IBS_z_0(p+1:mesh_num,p)   =0.5*kernel_mod_IBS_z(s_q,s(p),0);
            K_mat_IBS_z_1(p+1:mesh_num,p)   =0.5*kernel_mod_IBS_z(s_q,s(p),1);
            K_mat_IBS_z_0_pp(p+1:mesh_num,p)=0.5*kernel_mod_IBS_z(s_q,s(p),0.5);
            K_mat_IBS_z_1_pp(p+1:mesh_num,p)=0.5*kernel_mod_IBS_z(s_q,s(p),1.5);
            K_mat_IBS_z_2_pp(p+1:mesh_num,p)=0.5*kernel_mod_IBS_z(s_q,s(p),2.5);
            end
        else
            if (D_diffusion==1)
            K_mat_IBS_zz_1(p+1:mesh_num,p)  =kernel_mod_IBS_zz(s_q,s(p),1);
            K_mat_IBS_zz_2(p+1:mesh_num,p)  =kernel_mod_IBS_zz(s_q,s(p),2);
            K_mat_IBS_zz_30(p+1:mesh_num,p) =kernel_mod_IBS_zz(s_q,s(p),30);
            K_mat_IBS_zz_31(p+1:mesh_num,p) =kernel_mod_IBS_zz(s_q,s(p),31);
            K_mat_IBS_zz_4(p+1:mesh_num,p)  =kernel_mod_IBS_zz(s_q,s(p),4);
            end
            if (D_friction==1)
            K_mat_IBS_z_0(p+1:mesh_num,p)   =kernel_mod_IBS_z(s_q,s(p),0);
            K_mat_IBS_z_1(p+1:mesh_num,p)   =kernel_mod_IBS_z(s_q,s(p),1);
            K_mat_IBS_z_0_pp(p+1:mesh_num,p)=kernel_mod_IBS_z(s_q,s(p),0.5);
            K_mat_IBS_z_1_pp(p+1:mesh_num,p)=kernel_mod_IBS_z(s_q,s(p),1.5);
            K_mat_IBS_z_2_pp(p+1:mesh_num,p)=kernel_mod_IBS_z(s_q,s(p),2.5);
            end
        end
    end
    end
    %}
    % -------------- self-consistent (direct) solution ------------------ %
    
    if (iEnergy_gain_calc==0)
        if (iiterative==0)
            gkd_mat=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_k);
        else
            gkd_mat=iterative_mat_gen(K_mat,mesh_num,mesh_size,order_N)*transpose(G0_k);
        end
        G_c_11=gkd_mat/G0_k(1);
        G_11=abs(G_c_11);
        Gf_11(iscan)=G_11(end);
        max_G_11(iscan)=max(G_11);
        G_total(:,iscan)=G_11;
    elseif (iEnergy_gain_calc==1 && iTransverse_gain_calc==0)
        %{
        if (iiterative==0)
            one_minus_K_mat_inv=inv(eye(mesh_num)-mesh_size*K_mat); % using inv() sometimes does not work
        else
            one_minus_K_mat_inv=iterative_mat_gen(K_mat,mesh_num,mesh_size,order_N);
        end
        MBI_kernel_mat_11=one_minus_K_mat_inv;
        MBI_kernel_mat_12=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_13=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_14=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_21=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_22=one_minus_K_mat_inv;
        MBI_kernel_mat_23=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_24=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_31=mesh_size*M_minus_L_mat*one_minus_K_mat_inv;
        MBI_kernel_mat_32=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_33=eye(mesh_num);
        MBI_kernel_mat_34=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_41=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_42=MBI_kernel_mat_31;
        MBI_kernel_mat_43=zeros(mesh_num,mesh_num);
        MBI_kernel_mat_44=eye(mesh_num);
    
        MBI_kernel_mat=[MBI_kernel_mat_11 MBI_kernel_mat_12 MBI_kernel_mat_13 MBI_kernel_mat_14;...
                        MBI_kernel_mat_21 MBI_kernel_mat_22 MBI_kernel_mat_23 MBI_kernel_mat_24;...
                        MBI_kernel_mat_31 MBI_kernel_mat_32 MBI_kernel_mat_33 MBI_kernel_mat_34;...
                        MBI_kernel_mat_41 MBI_kernel_mat_42 MBI_kernel_mat_43 MBI_kernel_mat_44];
    
        tot_col_vec_old=[G0_k G0_kp E0_k E0_kp];
        tot_col_vec_old=transpose(tot_col_vec_old);
        %}
        %{
        % read tot_col_vec_old from external file
        prompt='read tot_col_vec_old from external files? (1:yes,0:no)';
        ans=input(prompt);
        if (ans==1)
            tot_col_vec_old=read_tot_col_vec(mesh_num,2);
        end
        %}
        %
        if (iIBS==0)
            tot_col_vec_new(1:mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_k);
            tot_col_vec_new(1*mesh_num+1:2*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_kp);
            tot_col_vec_new(2*mesh_num+1:3*mesh_num)=mesh_size*M_minus_L_mat*((eye(mesh_num)-mesh_size*K_mat)\transpose(G0_k))+transpose(E0_k);
            tot_col_vec_new(3*mesh_num+1:4*mesh_num)=mesh_size*M_minus_L_mat*((eye(mesh_num)-mesh_size*K_mat)\transpose(G0_kp))+transpose(E0_kp);
        elseif (iIBS==1)
            P_mat=eye(mesh_num)-mesh_size*K_mat-mesh_size*K_mat_IBS_z_1+2*mesh_size*K_mat_IBS_zz_2;
            Q_mat=-1i*mesh_size*K_mat_IBS_z_0_pp-1i*mesh_size*K_mat_IBS_zz_30;
            R_mat=-mesh_size*M_minus_L_mat-1i*mesh_size*K_mat_IBS_z_0-2*1i*mesh_size*K_mat_IBS_z_1_pp+4*1i*mesh_size*K_mat_IBS_zz_1-2*1i*K_mat_IBS_zz_31;
            S_mat=eye(mesh_num)+mesh_size*K_mat_IBS_z_0_pp-mesh_size*K_mat_IBS_z_2_pp+3*mesh_size*K_mat_IBS_zz_2-mesh_size*K_mat_IBS_zz_4;
            %{
            % below only consider diffusion coefficient based on Stupakov WEPSO68
            P_mat=eye(mesh_num)-mesh_size*K_mat+2*mesh_size*K_mat_IBS_zz_2;
            Q_mat=-1i*mesh_size*K_mat_IBS_zz_30;
            R_mat=-mesh_size*M_minus_L_mat+4*1i*mesh_size*K_mat_IBS_zz_1-2*1i*K_mat_IBS_zz_31;
            S_mat=eye(mesh_num)+3*mesh_size*K_mat_IBS_zz_2-mesh_size*K_mat_IBS_zz_4;
            %}
            %{
            % this does not work since det(Q_mat), det(R_mat) vanish
            Dpk0_vec=S_mat\transpose(E0_k)-Q_mat\(P_mat*(R_mat\transpose(E0_k)));
            Cbk0_vec=-S_mat\(R_mat*(P_mat\transpose(G0_k)))+Q_mat\transpose(G0_k);
            Bpk0_vec=-P_mat\(Q_mat*Dpk0_vec);
            Abk0_vec=P_mat\transpose(G0_k)-P_mat\(Q_mat*Cbk0_vec);
            %}
            %T_mat=S_mat-R_mat*inv(P_mat)*Q_mat;
            T_mat=S_mat-R_mat*(P_mat\Q_mat);
            Dpk0_vec=T_mat\transpose(E0_k);
            Cbk0_vec=-T_mat\(R_mat*(P_mat\transpose(G0_k)));
            Bpk0_vec=-P_mat\(Q_mat*(T_mat\transpose(E0_k)));
            Abk0_vec=P_mat\transpose(G0_k)-P_mat\(Q_mat*Cbk0_vec);
            
            tot_col_vec_new(1:mesh_num)=Abk0_vec+Bpk0_vec;
            tot_col_vec_new(1*mesh_num+1:2*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_kp);                                            % not including IBS yet
            tot_col_vec_new(2*mesh_num+1:3*mesh_num)=Cbk0_vec+Dpk0_vec;
            tot_col_vec_new(3*mesh_num+1:4*mesh_num)=mesh_size*M_minus_L_mat*((eye(mesh_num)-mesh_size*K_mat)\transpose(G0_kp))+transpose(E0_kp); % not including IBS yet
        end
        %}
        %tot_col_vec_new=MBI_kernel_mat*tot_col_vec_old;
        
    elseif (iEnergy_gain_calc==1 && iTransverse_gain_calc==1)
        if (iiterative==0)
            one_minus_K_mat_inv=inv(eye(mesh_num)-mesh_size*K_mat);
        else
            one_minus_K_mat_inv=iterative_mat_gen(K_mat,mesh_num,mesh_size,order_N);
        end
        
        % inv(I-K) matrix
        MBI_kernel_mat_01_01=one_minus_K_mat_inv; MBI_kernel_mat(0*mesh_num+1:1*mesh_num,0*mesh_num+1:1*mesh_num)=MBI_kernel_mat_01_01;
        MBI_kernel_mat_02_02=one_minus_K_mat_inv; MBI_kernel_mat(1*mesh_num+1:2*mesh_num,1*mesh_num+1:2*mesh_num)=MBI_kernel_mat_02_02;
        MBI_kernel_mat_03_03=one_minus_K_mat_inv; MBI_kernel_mat(2*mesh_num+1:3*mesh_num,2*mesh_num+1:3*mesh_num)=MBI_kernel_mat_03_03;
        MBI_kernel_mat_04_04=one_minus_K_mat_inv; MBI_kernel_mat(3*mesh_num+1:4*mesh_num,3*mesh_num+1:4*mesh_num)=MBI_kernel_mat_04_04;
        MBI_kernel_mat_05_05=one_minus_K_mat_inv; MBI_kernel_mat(4*mesh_num+1:5*mesh_num,4*mesh_num+1:5*mesh_num)=MBI_kernel_mat_05_05;
        MBI_kernel_mat_06_06=one_minus_K_mat_inv; MBI_kernel_mat(5*mesh_num+1:6*mesh_num,5*mesh_num+1:6*mesh_num)=MBI_kernel_mat_06_06;
        
        % (M-L)*inv(I-K) matrix
        MBI_kernel_mat_07_01=mesh_size*M_minus_L_mat*one_minus_K_mat_inv; MBI_kernel_mat(6*mesh_num+1:7*mesh_num,0*mesh_num+1:1*mesh_num)=MBI_kernel_mat_07_01;
        MBI_kernel_mat_08_02=mesh_size*M_minus_L_mat*one_minus_K_mat_inv; MBI_kernel_mat(7*mesh_num+1:8*mesh_num,1*mesh_num+1:2*mesh_num)=MBI_kernel_mat_08_02;
        MBI_kernel_mat_09_03=mesh_size*M_minus_L_mat*one_minus_K_mat_inv; MBI_kernel_mat(8*mesh_num+1:9*mesh_num,2*mesh_num+1:3*mesh_num)=MBI_kernel_mat_09_03;
        MBI_kernel_mat_10_04=mesh_size*M_minus_L_mat*one_minus_K_mat_inv; MBI_kernel_mat(9*mesh_num+1:10*mesh_num,3*mesh_num+1:4*mesh_num)=MBI_kernel_mat_10_04;
        MBI_kernel_mat_11_05=mesh_size*M_minus_L_mat*one_minus_K_mat_inv; MBI_kernel_mat(10*mesh_num+1:11*mesh_num,4*mesh_num+1:5*mesh_num)=MBI_kernel_mat_11_05;
        MBI_kernel_mat_12_06=mesh_size*M_minus_L_mat*one_minus_K_mat_inv; MBI_kernel_mat(11*mesh_num+1:12*mesh_num,5*mesh_num+1:6*mesh_num)=MBI_kernel_mat_12_06;
        
        % A*inv(I-K) matrix
        MBI_kernel_mat_13_01=mesh_size*A_mat*one_minus_K_mat_inv; MBI_kernel_mat(12*mesh_num+1:13*mesh_num,0*mesh_num+1:1*mesh_num)=MBI_kernel_mat_13_01;
        MBI_kernel_mat_14_02=mesh_size*A_mat*one_minus_K_mat_inv; MBI_kernel_mat(13*mesh_num+1:14*mesh_num,1*mesh_num+1:2*mesh_num)=MBI_kernel_mat_14_02;
        MBI_kernel_mat_15_03=mesh_size*A_mat*one_minus_K_mat_inv; MBI_kernel_mat(14*mesh_num+1:15*mesh_num,2*mesh_num+1:3*mesh_num)=MBI_kernel_mat_15_03;
        MBI_kernel_mat_16_04=mesh_size*A_mat*one_minus_K_mat_inv; MBI_kernel_mat(15*mesh_num+1:16*mesh_num,3*mesh_num+1:4*mesh_num)=MBI_kernel_mat_16_04;
        MBI_kernel_mat_17_05=mesh_size*A_mat*one_minus_K_mat_inv; MBI_kernel_mat(16*mesh_num+1:17*mesh_num,4*mesh_num+1:5*mesh_num)=MBI_kernel_mat_17_05;
        MBI_kernel_mat_18_06=mesh_size*A_mat*one_minus_K_mat_inv; MBI_kernel_mat(17*mesh_num+1:18*mesh_num,5*mesh_num+1:6*mesh_num)=MBI_kernel_mat_18_06;
        
        % B*inv(I-K) matrix
        MBI_kernel_mat_19_01=mesh_size*B_mat*one_minus_K_mat_inv; MBI_kernel_mat(18*mesh_num+1:19*mesh_num,0*mesh_num+1:1*mesh_num)=MBI_kernel_mat_19_01;
        MBI_kernel_mat_20_02=mesh_size*B_mat*one_minus_K_mat_inv; MBI_kernel_mat(19*mesh_num+1:20*mesh_num,1*mesh_num+1:2*mesh_num)=MBI_kernel_mat_20_02;
        MBI_kernel_mat_21_03=mesh_size*B_mat*one_minus_K_mat_inv; MBI_kernel_mat(20*mesh_num+1:21*mesh_num,2*mesh_num+1:3*mesh_num)=MBI_kernel_mat_21_03;
        MBI_kernel_mat_22_04=mesh_size*B_mat*one_minus_K_mat_inv; MBI_kernel_mat(21*mesh_num+1:22*mesh_num,3*mesh_num+1:4*mesh_num)=MBI_kernel_mat_22_04;
        MBI_kernel_mat_23_05=mesh_size*B_mat*one_minus_K_mat_inv; MBI_kernel_mat(22*mesh_num+1:23*mesh_num,4*mesh_num+1:5*mesh_num)=MBI_kernel_mat_23_05;
        MBI_kernel_mat_24_06=mesh_size*B_mat*one_minus_K_mat_inv; MBI_kernel_mat(23*mesh_num+1:24*mesh_num,5*mesh_num+1:6*mesh_num)=MBI_kernel_mat_24_06;
        
        % C*inv(I-K) matrix
        MBI_kernel_mat_25_01=mesh_size*C_mat*one_minus_K_mat_inv; MBI_kernel_mat(24*mesh_num+1:25*mesh_num,0*mesh_num+1:1*mesh_num)=MBI_kernel_mat_25_01;
        MBI_kernel_mat_26_02=mesh_size*C_mat*one_minus_K_mat_inv; MBI_kernel_mat(25*mesh_num+1:26*mesh_num,1*mesh_num+1:2*mesh_num)=MBI_kernel_mat_26_02;
        MBI_kernel_mat_27_03=mesh_size*C_mat*one_minus_K_mat_inv; MBI_kernel_mat(26*mesh_num+1:27*mesh_num,2*mesh_num+1:3*mesh_num)=MBI_kernel_mat_27_03;
        MBI_kernel_mat_28_04=mesh_size*C_mat*one_minus_K_mat_inv; MBI_kernel_mat(27*mesh_num+1:28*mesh_num,3*mesh_num+1:4*mesh_num)=MBI_kernel_mat_28_04;
        MBI_kernel_mat_29_05=mesh_size*C_mat*one_minus_K_mat_inv; MBI_kernel_mat(28*mesh_num+1:29*mesh_num,4*mesh_num+1:5*mesh_num)=MBI_kernel_mat_29_05;
        MBI_kernel_mat_30_06=mesh_size*C_mat*one_minus_K_mat_inv; MBI_kernel_mat(29*mesh_num+1:30*mesh_num,5*mesh_num+1:6*mesh_num)=MBI_kernel_mat_30_06;
        
        % D*inv(I-K) matrix
        MBI_kernel_mat_31_01=mesh_size*D_mat*one_minus_K_mat_inv; MBI_kernel_mat(30*mesh_num+1:31*mesh_num,0*mesh_num+1:1*mesh_num)=MBI_kernel_mat_31_01;
        MBI_kernel_mat_32_02=mesh_size*D_mat*one_minus_K_mat_inv; MBI_kernel_mat(31*mesh_num+1:32*mesh_num,1*mesh_num+1:2*mesh_num)=MBI_kernel_mat_32_02;
        MBI_kernel_mat_33_03=mesh_size*D_mat*one_minus_K_mat_inv; MBI_kernel_mat(32*mesh_num+1:33*mesh_num,2*mesh_num+1:3*mesh_num)=MBI_kernel_mat_33_03;
        MBI_kernel_mat_34_04=mesh_size*D_mat*one_minus_K_mat_inv; MBI_kernel_mat(33*mesh_num+1:34*mesh_num,3*mesh_num+1:4*mesh_num)=MBI_kernel_mat_34_04;
        MBI_kernel_mat_35_05=mesh_size*D_mat*one_minus_K_mat_inv; MBI_kernel_mat(34*mesh_num+1:35*mesh_num,4*mesh_num+1:5*mesh_num)=MBI_kernel_mat_35_05;
        MBI_kernel_mat_36_06=mesh_size*D_mat*one_minus_K_mat_inv; MBI_kernel_mat(35*mesh_num+1:36*mesh_num,5*mesh_num+1:6*mesh_num)=MBI_kernel_mat_36_06;
        
        clear MBI_kernel_mat;
        
        tot_col_vec_old=[G0_k G0_kp G0_kxz G0_kxpz G0_kyz G0_kypz...
                         E0_k E0_kp E0_kxz E0_kxpz E0_kyz E0_kypz...
                         Ax0_k Ax0_kp Ax0_kxz Ax0_kxpz Ax0_kyz Ax0_kypz...
                         Axp0_k Axp0_kp Axp0_kxz Axp0_kxpz Axp0_kyz Axp0_kypz...
                         Ay0_k Ay0_kp Ay0_kxz Ay0_kxpz Ay0_kyz Ay0_kypz...
                         Ayp0_k Ayp0_kp Ayp0_kxz Ayp0_kxpz Ayp0_kyz Ayp0_kypz];
        
        tot_col_vec_old=transpose(tot_col_vec_old);
        
        % direct use of the command below can run out of system memory, try
        % alternative method to multiply MBI_kernel_mat by tot_col_vec_old
        %tot_col_vec_new=MBI_kernel_mat*tot_col_vec_old;
        
        % use mtimesx, may speed but memory issue exists
        %tot_col_vec_new=mtimesx(MBI_kernel_mat,tot_col_vec_old);
        
        %{
        tot_col_vec_new(0*mesh_num+1:1*mesh_num)=mtimesx(one_minus_K_mat_inv,tot_col_vec_old(0*mesh_num+1:1*mesh_num));
        tot_col_vec_new(1*mesh_num+1:2*mesh_num)=mtimesx(one_minus_K_mat_inv,tot_col_vec_old(1*mesh_num+1:2*mesh_num));
        tot_col_vec_new(2*mesh_num+1:3*mesh_num)=mtimesx(one_minus_K_mat_inv,tot_col_vec_old(2*mesh_num+1:3*mesh_num));
        tot_col_vec_new(3*mesh_num+1:4*mesh_num)=mtimesx(one_minus_K_mat_inv,tot_col_vec_old(3*mesh_num+1:4*mesh_num));
        tot_col_vec_new(4*mesh_num+1:5*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(4*mesh_num+1:5*mesh_num));
        tot_col_vec_new(5*mesh_num+1:6*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(5*mesh_num+1:6*mesh_num));
        tot_col_vec_new(6*mesh_num+1:7*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(6*mesh_num+1:7*mesh_num))+mtimesx(MBI_kernel_mat_07_01,tot_col_vec_old(0*mesh_num+1:1*mesh_num));
        tot_col_vec_new(7*mesh_num+1:8*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(7*mesh_num+1:8*mesh_num))+mtimesx(MBI_kernel_mat_08_02,tot_col_vec_old(1*mesh_num+1:2*mesh_num));
        tot_col_vec_new(8*mesh_num+1:9*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(8*mesh_num+1:9*mesh_num))+mtimesx(MBI_kernel_mat_09_03,tot_col_vec_old(2*mesh_num+1:3*mesh_num));
        tot_col_vec_new(9*mesh_num+1:10*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(9*mesh_num+1:10*mesh_num))+mtimesx(MBI_kernel_mat_10_04,tot_col_vec_old(3*mesh_num+1:4*mesh_num));
        tot_col_vec_new(10*mesh_num+1:11*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(10*mesh_num+1:11*mesh_num))+mtimesx(MBI_kernel_mat_11_05,tot_col_vec_old(4*mesh_num+1:5*mesh_num));
        tot_col_vec_new(11*mesh_num+1:12*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(11*mesh_num+1:12*mesh_num))+mtimesx(MBI_kernel_mat_12_06,tot_col_vec_old(5*mesh_num+1:6*mesh_num));
        tot_col_vec_new(12*mesh_num+1:13*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(12*mesh_num+1:13*mesh_num))+mtimesx(MBI_kernel_mat_13_01,tot_col_vec_old(0*mesh_num+1:1*mesh_num));
        tot_col_vec_new(13*mesh_num+1:14*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(13*mesh_num+1:14*mesh_num))+mtimesx(MBI_kernel_mat_14_02,tot_col_vec_old(1*mesh_num+1:2*mesh_num));
        tot_col_vec_new(14*mesh_num+1:15*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(14*mesh_num+1:15*mesh_num))+mtimesx(MBI_kernel_mat_15_03,tot_col_vec_old(2*mesh_num+1:3*mesh_num));
        tot_col_vec_new(15*mesh_num+1:16*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(15*mesh_num+1:16*mesh_num))+mtimesx(MBI_kernel_mat_16_04,tot_col_vec_old(3*mesh_num+1:4*mesh_num));
        tot_col_vec_new(16*mesh_num+1:17*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(16*mesh_num+1:17*mesh_num))+mtimesx(MBI_kernel_mat_17_05,tot_col_vec_old(4*mesh_num+1:5*mesh_num));
        tot_col_vec_new(17*mesh_num+1:18*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(17*mesh_num+1:18*mesh_num))+mtimesx(MBI_kernel_mat_18_06,tot_col_vec_old(5*mesh_num+1:6*mesh_num));
        tot_col_vec_new(18*mesh_num+1:19*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(18*mesh_num+1:19*mesh_num))+mtimesx(MBI_kernel_mat_19_01,tot_col_vec_old(0*mesh_num+1:1*mesh_num));
        tot_col_vec_new(19*mesh_num+1:20*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(19*mesh_num+1:20*mesh_num))+mtimesx(MBI_kernel_mat_20_02,tot_col_vec_old(1*mesh_num+1:2*mesh_num));
        tot_col_vec_new(20*mesh_num+1:21*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(20*mesh_num+1:21*mesh_num))+mtimesx(MBI_kernel_mat_21_03,tot_col_vec_old(2*mesh_num+1:3*mesh_num));
        tot_col_vec_new(21*mesh_num+1:22*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(21*mesh_num+1:22*mesh_num))+mtimesx(MBI_kernel_mat_22_04,tot_col_vec_old(3*mesh_num+1:4*mesh_num));
        tot_col_vec_new(22*mesh_num+1:23*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(22*mesh_num+1:23*mesh_num))+mtimesx(MBI_kernel_mat_23_05,tot_col_vec_old(4*mesh_num+1:5*mesh_num));
        tot_col_vec_new(23*mesh_num+1:24*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(23*mesh_num+1:24*mesh_num))+mtimesx(MBI_kernel_mat_24_06,tot_col_vec_old(5*mesh_num+1:6*mesh_num));
        tot_col_vec_new(24*mesh_num+1:25*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(24*mesh_num+1:25*mesh_num))+mtimesx(MBI_kernel_mat_25_01,tot_col_vec_old(0*mesh_num+1:1*mesh_num));
        tot_col_vec_new(25*mesh_num+1:26*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(25*mesh_num+1:26*mesh_num))+mtimesx(MBI_kernel_mat_26_02,tot_col_vec_old(1*mesh_num+1:2*mesh_num));
        tot_col_vec_new(26*mesh_num+1:27*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(26*mesh_num+1:27*mesh_num))+mtimesx(MBI_kernel_mat_27_03,tot_col_vec_old(2*mesh_num+1:3*mesh_num));
        tot_col_vec_new(27*mesh_num+1:28*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(27*mesh_num+1:28*mesh_num))+mtimesx(MBI_kernel_mat_28_04,tot_col_vec_old(3*mesh_num+1:4*mesh_num));
        tot_col_vec_new(28*mesh_num+1:29*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(28*mesh_num+1:29*mesh_num))+mtimesx(MBI_kernel_mat_29_05,tot_col_vec_old(4*mesh_num+1:5*mesh_num));
        tot_col_vec_new(29*mesh_num+1:30*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(29*mesh_num+1:30*mesh_num))+mtimesx(MBI_kernel_mat_30_06,tot_col_vec_old(5*mesh_num+1:6*mesh_num));
        tot_col_vec_new(30*mesh_num+1:31*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(30*mesh_num+1:31*mesh_num))+mtimesx(MBI_kernel_mat_31_01,tot_col_vec_old(0*mesh_num+1:1*mesh_num));
        tot_col_vec_new(31*mesh_num+1:32*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(31*mesh_num+1:32*mesh_num))+mtimesx(MBI_kernel_mat_32_02,tot_col_vec_old(1*mesh_num+1:2*mesh_num));
        tot_col_vec_new(32*mesh_num+1:33*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(32*mesh_num+1:33*mesh_num))+mtimesx(MBI_kernel_mat_33_03,tot_col_vec_old(2*mesh_num+1:3*mesh_num));
        tot_col_vec_new(33*mesh_num+1:34*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(33*mesh_num+1:34*mesh_num))+mtimesx(MBI_kernel_mat_34_04,tot_col_vec_old(3*mesh_num+1:4*mesh_num));
        tot_col_vec_new(34*mesh_num+1:35*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(34*mesh_num+1:35*mesh_num))+mtimesx(MBI_kernel_mat_35_05,tot_col_vec_old(4*mesh_num+1:5*mesh_num));
        tot_col_vec_new(35*mesh_num+1:36*mesh_num)=mtimesx(eye(mesh_num),tot_col_vec_old(35*mesh_num+1:36*mesh_num))+mtimesx(MBI_kernel_mat_36_06,tot_col_vec_old(5*mesh_num+1:6*mesh_num));
        %}
        %{
        % 20200704, work in progress
        tot_col_vec_new(0*mesh_num+1:1*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(1*mesh_num+1:2*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(2*mesh_num+1:3*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(3*mesh_num+1:4*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(4*mesh_num+1:5*mesh_num)=eye(mesh_num)*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(5*mesh_num+1:6*mesh_num)=eye(mesh_num)*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        %}
        %
        tot_col_vec_new(0*mesh_num+1:1*mesh_num)=one_minus_K_mat_inv*tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(1*mesh_num+1:2*mesh_num)=one_minus_K_mat_inv*tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(2*mesh_num+1:3*mesh_num)=one_minus_K_mat_inv*tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(3*mesh_num+1:4*mesh_num)=one_minus_K_mat_inv*tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(4*mesh_num+1:5*mesh_num)=eye(mesh_num)*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(5*mesh_num+1:6*mesh_num)=eye(mesh_num)*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        tot_col_vec_new(6*mesh_num+1:7*mesh_num)=eye(mesh_num)*tot_col_vec_old(6*mesh_num+1:7*mesh_num)+MBI_kernel_mat_07_01*tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(7*mesh_num+1:8*mesh_num)=eye(mesh_num)*tot_col_vec_old(7*mesh_num+1:8*mesh_num)+MBI_kernel_mat_08_02*tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(8*mesh_num+1:9*mesh_num)=eye(mesh_num)*tot_col_vec_old(8*mesh_num+1:9*mesh_num)+MBI_kernel_mat_09_03*tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(9*mesh_num+1:10*mesh_num)=eye(mesh_num)*tot_col_vec_old(9*mesh_num+1:10*mesh_num)+MBI_kernel_mat_10_04*tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(10*mesh_num+1:11*mesh_num)=eye(mesh_num)*tot_col_vec_old(10*mesh_num+1:11*mesh_num)+MBI_kernel_mat_11_05*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(11*mesh_num+1:12*mesh_num)=eye(mesh_num)*tot_col_vec_old(11*mesh_num+1:12*mesh_num)+MBI_kernel_mat_12_06*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        tot_col_vec_new(12*mesh_num+1:13*mesh_num)=eye(mesh_num)*tot_col_vec_old(12*mesh_num+1:13*mesh_num)+MBI_kernel_mat_13_01*tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(13*mesh_num+1:14*mesh_num)=eye(mesh_num)*tot_col_vec_old(13*mesh_num+1:14*mesh_num)+MBI_kernel_mat_14_02*tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(14*mesh_num+1:15*mesh_num)=eye(mesh_num)*tot_col_vec_old(14*mesh_num+1:15*mesh_num)+MBI_kernel_mat_15_03*tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(15*mesh_num+1:16*mesh_num)=eye(mesh_num)*tot_col_vec_old(15*mesh_num+1:16*mesh_num)+MBI_kernel_mat_16_04*tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(16*mesh_num+1:17*mesh_num)=eye(mesh_num)*tot_col_vec_old(16*mesh_num+1:17*mesh_num)+MBI_kernel_mat_17_05*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(17*mesh_num+1:18*mesh_num)=eye(mesh_num)*tot_col_vec_old(17*mesh_num+1:18*mesh_num)+MBI_kernel_mat_18_06*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        tot_col_vec_new(18*mesh_num+1:19*mesh_num)=eye(mesh_num)*tot_col_vec_old(18*mesh_num+1:19*mesh_num)+MBI_kernel_mat_19_01*tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(19*mesh_num+1:20*mesh_num)=eye(mesh_num)*tot_col_vec_old(19*mesh_num+1:20*mesh_num)+MBI_kernel_mat_20_02*tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(20*mesh_num+1:21*mesh_num)=eye(mesh_num)*tot_col_vec_old(20*mesh_num+1:21*mesh_num)+MBI_kernel_mat_21_03*tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(21*mesh_num+1:22*mesh_num)=eye(mesh_num)*tot_col_vec_old(21*mesh_num+1:22*mesh_num)+MBI_kernel_mat_22_04*tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(22*mesh_num+1:23*mesh_num)=eye(mesh_num)*tot_col_vec_old(22*mesh_num+1:23*mesh_num)+MBI_kernel_mat_23_05*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(23*mesh_num+1:24*mesh_num)=eye(mesh_num)*tot_col_vec_old(23*mesh_num+1:24*mesh_num)+MBI_kernel_mat_24_06*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        tot_col_vec_new(24*mesh_num+1:25*mesh_num)=eye(mesh_num)*tot_col_vec_old(24*mesh_num+1:25*mesh_num)+MBI_kernel_mat_25_01*tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(25*mesh_num+1:26*mesh_num)=eye(mesh_num)*tot_col_vec_old(25*mesh_num+1:26*mesh_num)+MBI_kernel_mat_26_02*tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(26*mesh_num+1:27*mesh_num)=eye(mesh_num)*tot_col_vec_old(26*mesh_num+1:27*mesh_num)+MBI_kernel_mat_27_03*tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(27*mesh_num+1:28*mesh_num)=eye(mesh_num)*tot_col_vec_old(27*mesh_num+1:28*mesh_num)+MBI_kernel_mat_28_04*tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(28*mesh_num+1:29*mesh_num)=eye(mesh_num)*tot_col_vec_old(28*mesh_num+1:29*mesh_num)+MBI_kernel_mat_29_05*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(29*mesh_num+1:30*mesh_num)=eye(mesh_num)*tot_col_vec_old(29*mesh_num+1:30*mesh_num)+MBI_kernel_mat_30_06*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        tot_col_vec_new(30*mesh_num+1:31*mesh_num)=eye(mesh_num)*tot_col_vec_old(30*mesh_num+1:31*mesh_num)+MBI_kernel_mat_31_01*tot_col_vec_old(0*mesh_num+1:1*mesh_num);
        tot_col_vec_new(31*mesh_num+1:32*mesh_num)=eye(mesh_num)*tot_col_vec_old(31*mesh_num+1:32*mesh_num)+MBI_kernel_mat_32_02*tot_col_vec_old(1*mesh_num+1:2*mesh_num);
        tot_col_vec_new(32*mesh_num+1:33*mesh_num)=eye(mesh_num)*tot_col_vec_old(32*mesh_num+1:33*mesh_num)+MBI_kernel_mat_33_03*tot_col_vec_old(2*mesh_num+1:3*mesh_num);
        tot_col_vec_new(33*mesh_num+1:34*mesh_num)=eye(mesh_num)*tot_col_vec_old(33*mesh_num+1:34*mesh_num)+MBI_kernel_mat_34_04*tot_col_vec_old(3*mesh_num+1:4*mesh_num);
        tot_col_vec_new(34*mesh_num+1:35*mesh_num)=eye(mesh_num)*tot_col_vec_old(34*mesh_num+1:35*mesh_num)+MBI_kernel_mat_35_05*tot_col_vec_old(4*mesh_num+1:5*mesh_num);
        tot_col_vec_new(35*mesh_num+1:36*mesh_num)=eye(mesh_num)*tot_col_vec_old(35*mesh_num+1:36*mesh_num)+MBI_kernel_mat_36_06*tot_col_vec_old(5*mesh_num+1:6*mesh_num);
        %}
        
        % clear large-size variables; release system memory
        clear MBI_kernel_mat_01_01 MBI_kernel_mat_02_02 MBI_kernel_mat_03_03 MBI_kernel_mat_04_04 MBI_kernel_mat_05_05 MBI_kernel_mat_06_06
        clear MBI_kernel_mat_07_01 MBI_kernel_mat_08_02 MBI_kernel_mat_09_03 MBI_kernel_mat_10_04 MBI_kernel_mat_11_05 MBI_kernel_mat_12_06
        clear MBI_kernel_mat_13_01 MBI_kernel_mat_14_02 MBI_kernel_mat_15_03 MBI_kernel_mat_16_04 MBI_kernel_mat_17_05 MBI_kernel_mat_18_06
        clear MBI_kernel_mat_19_01 MBI_kernel_mat_20_02 MBI_kernel_mat_21_03 MBI_kernel_mat_22_04 MBI_kernel_mat_23_05 MBI_kernel_mat_24_06
        clear MBI_kernel_mat_25_01 MBI_kernel_mat_26_02 MBI_kernel_mat_27_03 MBI_kernel_mat_28_04 MBI_kernel_mat_29_05 MBI_kernel_mat_30_06
        clear MBI_kernel_mat_31_01 MBI_kernel_mat_32_02 MBI_kernel_mat_33_03 MBI_kernel_mat_34_04 MBI_kernel_mat_35_05 MBI_kernel_mat_36_06
        
    end
    
    if (iEnergy_gain_calc==1)       
        if (iTransverse_gain_calc==0)
            gkd_mat=tot_col_vec_new(1:mesh_num);
            gkp_mat=tot_col_vec_new(mesh_num+1:2*mesh_num); 
            ekd_mat=tot_col_vec_new(2*mesh_num+1:3*mesh_num);
            ekp_mat=tot_col_vec_new(3*mesh_num+1:4*mesh_num);
        end
        
        if (iTransverse_gain_calc==1)
            gkd_mat=tot_col_vec_new(0*mesh_num+1:1*mesh_num);
            gkp_mat=tot_col_vec_new(1*mesh_num+1:2*mesh_num); 
            gkxz_mat=tot_col_vec_new(2*mesh_num+1:3*mesh_num);
            gkxpz_mat=tot_col_vec_new(3*mesh_num+1:4*mesh_num);
            gkyz_mat=tot_col_vec_new(4*mesh_num+1:5*mesh_num);
            gkypz_mat=tot_col_vec_new(5*mesh_num+1:6*mesh_num);
            
            ekd_mat=tot_col_vec_new(6*mesh_num+1:7*mesh_num);
            ekp_mat=tot_col_vec_new(7*mesh_num+1:8*mesh_num);
            ekxz_mat=tot_col_vec_new(8*mesh_num+1:9*mesh_num);
            ekxpz_mat=tot_col_vec_new(9*mesh_num+1:10*mesh_num);
            ekyz_mat=tot_col_vec_new(10*mesh_num+1:11*mesh_num);
            ekypz_mat=tot_col_vec_new(11*mesh_num+1:12*mesh_num);
            
            axkd_mat=tot_col_vec_new(12*mesh_num+1:13*mesh_num);
            axkp_mat=tot_col_vec_new(13*mesh_num+1:14*mesh_num);
            axkxz_mat=tot_col_vec_new(14*mesh_num+1:15*mesh_num);
            axkxpz_mat=tot_col_vec_new(15*mesh_num+1:16*mesh_num);
            axkyz_mat=tot_col_vec_new(16*mesh_num+1:17*mesh_num);
            axkypz_mat=tot_col_vec_new(17*mesh_num+1:18*mesh_num);
            
            axpkd_mat=tot_col_vec_new(18*mesh_num+1:19*mesh_num);
            axpkp_mat=tot_col_vec_new(19*mesh_num+1:20*mesh_num);
            axpkxz_mat=tot_col_vec_new(20*mesh_num+1:21*mesh_num);
            axpkxpz_mat=tot_col_vec_new(21*mesh_num+1:22*mesh_num);
            axpkyz_mat=tot_col_vec_new(22*mesh_num+1:23*mesh_num);
            axpkypz_mat=tot_col_vec_new(23*mesh_num+1:24*mesh_num);
            
            aykd_mat=tot_col_vec_new(24*mesh_num+1:25*mesh_num);
            aykp_mat=tot_col_vec_new(25*mesh_num+1:26*mesh_num);
            aykxz_mat=tot_col_vec_new(26*mesh_num+1:27*mesh_num);
            aykxpz_mat=tot_col_vec_new(27*mesh_num+1:28*mesh_num);
            aykyz_mat=tot_col_vec_new(28*mesh_num+1:29*mesh_num);
            aykypz_mat=tot_col_vec_new(29*mesh_num+1:30*mesh_num);
            
            aypkd_mat=tot_col_vec_new(30*mesh_num+1:31*mesh_num);
            aypkp_mat=tot_col_vec_new(31*mesh_num+1:32*mesh_num);
            aypkxz_mat=tot_col_vec_new(32*mesh_num+1:33*mesh_num);
            aypkxpz_mat=tot_col_vec_new(33*mesh_num+1:34*mesh_num);
            aypkyz_mat=tot_col_vec_new(34*mesh_num+1:35*mesh_num);
            aypkypz_mat=tot_col_vec_new(35*mesh_num+1:36*mesh_num);
        end
        
        G_c_11=gkd_mat/G0_k(1); G_11=abs(G_c_11); Gf_11(iscan)=G_11(end); max_G_11(iscan)=max(G_11);
        %G_c_22=gkp_mat/G0_kp(1);G_22=abs(G_c_22); Gf_22(iscan)=G_22(end);
        G_c_22=gkp_mat/e_1k0;        G_22=abs(G_c_22); Gf_22(iscan)=G_22(end);
        G_c_31=ekd_mat/G0_k(1); G_31=abs(G_c_31); Gf_31(iscan)=G_31(end);
        %G_c_33=ekd_mat/E0_k(1); G_33=abs(G_c_33); Gf_33(iscan)=G_33(end);
        %G_c_42=ekp_mat/G0_kp(1);G_42=abs(G_c_42); Gf_42(iscan)=G_42(end);
        %G_c_44=ekp_mat/E0_kp(1);G_44=abs(G_c_44); Gf_44(iscan)=G_44(end);
        G_c_44=ekp_mat/e_1k0;        G_44=abs(G_c_44); Gf_44(iscan)=G_44(end);
        
        G_total(:,iscan)=abs(G_c_11+G_c_22);
        Gc_total(:,iscan)=G_c_11;               % 20200920 for test of SES_calc
        E_total(:,iscan)=abs(G_c_31+G_c_44); 
        Gf_total(iscan)=abs(G_c_11(end)+G_c_22(end));
        Ef_total(iscan)=abs(G_c_31(end)+G_c_44(end));
        
        gkd_mat_fin_abs(iscan)=abs(gkd_mat(end));
        gkp_mat_fin_abs(iscan)=abs(gkp_mat(end));
        ekd_mat_fin_abs(iscan)=abs(ekd_mat(end));
        ekp_mat_fin_abs(iscan)=abs(ekp_mat(end));
        
        gkd_mat_fin(iscan)=gkd_mat(end);
        gkp_mat_fin(iscan)=gkp_mat(end);
        ekd_mat_fin(iscan)=ekd_mat(end);
        ekp_mat_fin(iscan)=ekp_mat(end);
        
        gk_fin(iscan)=gkd_mat(end)+gkp_mat(end);
        ek_fin(iscan)=ekd_mat(end)+ekp_mat(end);
        
        if (iTransverse_gain_calc==1)
            % absolute value
            gkxz_mat_fin_abs(iscan)=abs(gkxz_mat(end));
            gkxpz_mat_fin_abs(iscan)=abs(gkxpz_mat(end));
            gkyz_mat_fin_abs(iscan)=abs(gkyz_mat(end));
            gkypz_mat_fin_abs(iscan)=abs(gkypz_mat(end));
            ekxz_mat_fin_abs(iscan)=abs(ekxz_mat(end));
            ekxpz_mat_fin_abs(iscan)=abs(ekxpz_mat(end));
            ekyz_mat_fin_abs(iscan)=abs(ekyz_mat(end));
            ekypz_mat_fin_abs(iscan)=abs(ekypz_mat(end));
            axkd_mat_fin_abs(iscan)=abs(axkd_mat(end));
            axkp_mat_fin_abs(iscan)=abs(axkp_mat(end));
            axkxz_mat_fin_abs(iscan)=abs(axkxz_mat(end));
            axkxpz_mat_fin_abs(iscan)=abs(axkxpz_mat(end));
            axkyz_mat_fin_abs(iscan)=abs(axkyz_mat(end));
            axkypz_mat_fin_abs(iscan)=abs(axkypz_mat(end));
            axpkd_mat_fin_abs(iscan)=abs(axpkd_mat(end));
            axpkp_mat_fin_abs(iscan)=abs(axpkp_mat(end));
            axpkxz_mat_fin_abs(iscan)=abs(axpkxz_mat(end));
            axpkxpz_mat_fin_abs(iscan)=abs(axpkxpz_mat(end));
            axpkyz_mat_fin_abs(iscan)=abs(axpkyz_mat(end));
            axpkypz_mat_fin_abs(iscan)=abs(axpkypz_mat(end));
            aykd_mat_fin_abs(iscan)=abs(aykd_mat(end));
            aykp_mat_fin_abs(iscan)=abs(aykp_mat(end));
            aykxz_mat_fin_abs(iscan)=abs(aykxz_mat(end));
            aykxpz_mat_fin_abs(iscan)=abs(aykxpz_mat(end));
            aykyz_mat_fin_abs(iscan)=abs(aykyz_mat(end));
            aykypz_mat_fin_abs(iscan)=abs(aykypz_mat(end));
            aypkd_mat_fin_abs(iscan)=abs(aypkd_mat(end));
            aypkp_mat_fin_abs(iscan)=abs(aypkp_mat(end));
            aypkxz_mat_fin_abs(iscan)=abs(aypkxz_mat(end));
            aypkxpz_mat_fin_abs(iscan)=abs(aypkxpz_mat(end));
            aypkyz_mat_fin_abs(iscan)=abs(aypkyz_mat(end));
            aypkypz_mat_fin_abs(iscan)=abs(aypkypz_mat(end));
            
            % complex value
            gkxz_mat_fin(iscan)=gkxz_mat(end);
            gkxpz_mat_fin(iscan)=gkxpz_mat(end);
            gkyz_mat_fin(iscan)=gkyz_mat(end);
            gkypz_mat_fin(iscan)=gkypz_mat(end);
            ekxz_mat_fin(iscan)=ekxz_mat(end);
            ekxpz_mat_fin(iscan)=ekxpz_mat(end);
            ekyz_mat_fin(iscan)=ekyz_mat(end);
            ekypz_mat_fin(iscan)=ekypz_mat(end);
            axkd_mat_fin(iscan)=axkd_mat(end);
            axkp_mat_fin(iscan)=axkp_mat(end);
            axkxz_mat_fin(iscan)=axkxz_mat(end);
            axkxpz_mat_fin(iscan)=axkxpz_mat(end);
            axkyz_mat_fin(iscan)=axkyz_mat(end);
            axkypz_mat_fin(iscan)=axkypz_mat(end);
            axpkd_mat_fin(iscan)=axpkd_mat(end);
            axpkp_mat_fin(iscan)=axpkp_mat(end);
            axpkxz_mat_fin(iscan)=axpkxz_mat(end);
            axpkxpz_mat_fin(iscan)=axpkxpz_mat(end);
            axpkyz_mat_fin(iscan)=axpkyz_mat(end);
            axpkypz_mat_fin(iscan)=axpkypz_mat(end);
            aykd_mat_fin(iscan)=aykd_mat(end);
            aykp_mat_fin(iscan)=aykp_mat(end);
            aykxz_mat_fin(iscan)=aykxz_mat(end);
            aykxpz_mat_fin(iscan)=aykxpz_mat(end);
            aykyz_mat_fin(iscan)=aykyz_mat(end);
            aykypz_mat_fin(iscan)=aykypz_mat(end);
            aypkd_mat_fin(iscan)=aypkd_mat(end);
            aypkp_mat_fin(iscan)=aypkp_mat(end);
            aypkxz_mat_fin(iscan)=aypkxz_mat(end);
            aypkxpz_mat_fin(iscan)=aypkxpz_mat(end);
            aypkyz_mat_fin(iscan)=aypkyz_mat(end);
            aypkypz_mat_fin(iscan)=aypkypz_mat(end);
            
            % physical values
            gk_fin(iscan)=gkd_mat(end)+gkp_mat(end)+gkxz_mat(end)+gkxpz_mat(end)+gkyz_mat(end)+gkypz_mat(end);
            ek_fin(iscan)=ekd_mat(end)+ekp_mat(end)+ekxz_mat(end)+ekxpz_mat(end)+ekyz_mat(end)+ekypz_mat(end);
            axk_fin(iscan)=axkd_mat(end)+axkp_mat(end)+axkxz_mat(end)+axkxpz_mat(end)+axkyz_mat(end)+axkypz_mat(end);
            axpk_fin(iscan)=axpkd_mat(end)+axpkp_mat(end)+axpkxz_mat(end)+axpkxpz_mat(end)+axpkyz_mat(end)+axpkypz_mat(end);
            ayk_fin(iscan)=aykd_mat(end)+aykp_mat(end)+aykxz_mat(end)+aykxpz_mat(end)+aykyz_mat(end)+aykypz_mat(end);
            aypk_fin(iscan)=aypkd_mat(end)+aypkp_mat(end)+aypkxz_mat(end)+aypkxpz_mat(end)+aypkyz_mat(end)+aypkypz_mat(end);
        end
        
        if (first_harmonic_notification==1)
            %if (iscan==1)
            if (iTransverse_gain_calc==1)
                LHS_01(:,iscan)=abs(interp1(s_ele,R56_ele,s)).*abs(ekd_mat);          % unit: cm
                LHS_02(:,iscan)=abs(interp1(s_ele,R56_ele,s)).*abs(ekp_mat);          % unit: cm
                LHS_03(:,iscan)=abs(interp1(s_ele,R56_ele,s)).*abs(ekd_mat+ekp_mat);  % unit: cm
                RHS(:,iscan)=lambda./interp1(s_ele,C_ele,s);                          % unit: cm
            else
                LHS_01(:,iscan)=abs(interp1(s_ele,R56_ele,s))'.*abs(ekd_mat);         % unit: cm
                LHS_02(:,iscan)=abs(interp1(s_ele,R56_ele,s))'.*abs(ekp_mat);         % unit: cm
                LHS_03(:,iscan)=abs(interp1(s_ele,R56_ele,s))'.*abs(ekd_mat+ekp_mat); % unit: cm
                RHS(:,iscan)=lambda./interp1(s_ele,C_ele,s);                          % unit: cm
            end
            %end
        end
        
    end
    
    if (iEnergy_gain_calc==1)
    	fprintf('iscan = %d, lambda = %.2f um, enx = %.2f um, eny = %.2f um, dp/p = %.2f, d-d gain = %f, d-e gain = %f, e-d gain = %f, e-e gain = %f \n',iscan,lambda*10^4,emit_norm_x*10^4,emit_norm_y*10^4,sigma_delta,Gf_11(iscan),Gf_31(iscan),Gf_22(iscan),Gf_44(iscan));
        if (iTransverse_gain_calc==1)
            fprintf('iscan = %d, lambda = %.2f um, enx = %.2f um, eny = %.2f um, dp/p = %.2f, bk_f = %f, pk_f = %f, axk_f = %f, axpk_f = %f, ayk_f = %f, aypk_f = %f \n',iscan,lambda*10^4,emit_norm_x*10^4,emit_norm_y*10^4,sigma_delta,abs(gk_fin(iscan)),abs(ek_fin(iscan)),abs(axk_fin(iscan)),abs(axpk_fin(iscan)),abs(ayk_fin(iscan)),abs(aypk_fin(iscan)));
        end
    else
    	fprintf('iscan = %d, lambda = %.2f um, enx = %.2f um, eny = %.2f um, dp/p = %.2f, density gain = %f \n',iscan,lambda*10^4,emit_norm_x*10^4,emit_norm_y*10^4,sigma_delta,Gf_11(iscan));
    end
    %{
    if (scan_num>=1 && iEnergy_gain_calc==1 && iTransverse_gain_calc==0)   
        if (ans_tmp==1) % write data to external file for next run
            Re_gkd=real(gkd_mat(end)); Im_gkd=imag(gkd_mat(end)); gkd_mat_mat=[Re_gkd Im_gkd];
            Re_gkp=real(gkp_mat(end)); Im_gkp=imag(gkp_mat(end)); gkp_mat_mat=[Re_gkp Im_gkp];
            Re_ekd=real(ekd_mat(end)); Im_ekd=imag(ekd_mat(end)); ekd_mat_mat=[Re_ekd Im_ekd];
            Re_ekp=real(ekp_mat(end)); Im_ekp=imag(ekp_mat(end)); ekp_mat_mat=[Re_ekp Im_ekp];
        
            if (iscan==1)
                file_gkd=fopen(sprintf('gkd_mat_%02d.bin',file_index_write),'a');
                file_gkp=fopen(sprintf('gkp_mat_%02d.bin',file_index_write),'a');
                file_ekd=fopen(sprintf('ekd_mat_%02d.bin',file_index_write),'a');
                file_ekp=fopen(sprintf('ekp_mat_%02d.bin',file_index_write),'a');
            end
            
            % format:
            % [Re(1) Im(1) Re(2) Im(2) Re(3) Im(3)...]
            fwrite(file_gkd,gkd_mat_mat,'double');
            fwrite(file_gkp,gkp_mat_mat,'double');
            fwrite(file_ekd,ekd_mat_mat,'double');
            fwrite(file_ekp,ekp_mat_mat,'double');
        
            if (iscan==scan_num)
                fclose(file_gkd);
                fclose(file_gkp);
                fclose(file_ekd);
                fclose(file_ekp);
            end
        end
    elseif (scan_num>=1 && iEnergy_gain_calc==1 && iTransverse_gain_calc==1)
        if (ans_tmp==1) % write data to external file for next run
            Re_gkd=real(gkd_mat(end)); Im_gkd=imag(gkd_mat(end)); gkd_mat_mat=[Re_gkd Im_gkd];
            Re_gkp=real(gkp_mat(end)); Im_gkp=imag(gkp_mat(end)); gkp_mat_mat=[Re_gkp Im_gkp];
            Re_gkxz=real(gkxz_mat(end)); Im_gkxz=imag(gkxz_mat(end)); gkxz_mat_mat=[Re_gkxz Im_gkxz];
            Re_gkxpz=real(gkxpz_mat(end)); Im_gkxpz=imag(gkxpz_mat(end)); gkxpz_mat_mat=[Re_gkxpz Im_gkxpz];
            Re_gkyz=real(gkyz_mat(end)); Im_gkyz=imag(gkyz_mat(end)); gkyz_mat_mat=[Re_gkyz Im_gkyz];
            Re_gkypz=real(gkypz_mat(end)); Im_gkypz=imag(gkypz_mat(end)); gkypz_mat_mat=[Re_gkypz Im_gkypz];
            
            Re_ekd=real(ekd_mat(end)); Im_ekd=imag(ekd_mat(end)); ekd_mat_mat=[Re_ekd Im_ekd];
            Re_ekp=real(ekp_mat(end)); Im_ekp=imag(ekp_mat(end)); ekp_mat_mat=[Re_ekp Im_ekp];            
            Re_ekxz=real(ekxz_mat(end)); Im_ekxz=imag(ekxz_mat(end)); ekxz_mat_mat=[Re_ekxz Im_ekxz];
            Re_ekxpz=real(ekxpz_mat(end)); Im_ekxpz=imag(ekxpz_mat(end)); ekxpz_mat_mat=[Re_ekxpz Im_ekxpz];
            Re_ekyz=real(ekyz_mat(end)); Im_ekyz=imag(ekyz_mat(end)); ekyz_mat_mat=[Re_ekyz Im_ekyz];
            Re_ekypz=real(ekypz_mat(end)); Im_ekypz=imag(ekypz_mat(end)); ekypz_mat_mat=[Re_ekypz Im_ekypz];
            
            Re_axkd=real(axkd_mat(end)); Im_axkd=imag(axkd_mat(end)); axkd_mat_mat=[Re_axkd Im_axkd];
            Re_axkp=real(axkp_mat(end)); Im_axkp=imag(axkp_mat(end)); axkp_mat_mat=[Re_axkp Im_axkp];
            Re_axkxz=real(axkxz_mat(end)); Im_axkxz=imag(axkxz_mat(end)); axkxz_mat_mat=[Re_axkxz Im_axkxz];
            Re_axkxpz=real(axkxpz_mat(end)); Im_axkxpz=imag(axkxpz_mat(end)); axkxpz_mat_mat=[Re_axkxpz Im_axkxpz];
            Re_axkyz=real(axkyz_mat(end)); Im_axkyz=imag(axkyz_mat(end)); axkyz_mat_mat=[Re_axkyz Im_axkyz];
            Re_axkypz=real(axkypz_mat(end)); Im_axkypz=imag(axkypz_mat(end)); axkypz_mat_mat=[Re_axkypz Im_axkypz];
            
            Re_axpkd=real(axpkd_mat(end)); Im_axpkd=imag(axpkd_mat(end)); axpkd_mat_mat=[Re_axpkd Im_axpkd];
            Re_axpkp=real(axpkp_mat(end)); Im_axpkp=imag(axpkp_mat(end)); axpkp_mat_mat=[Re_axpkp Im_axpkp];    
            Re_axpkxz=real(axpkxz_mat(end)); Im_axpkxz=imag(axpkxz_mat(end)); axpkxz_mat_mat=[Re_axpkxz Im_axpkxz];
            Re_axpkxpz=real(axpkxpz_mat(end)); Im_axpkxpz=imag(axpkxpz_mat(end)); axpkxpz_mat_mat=[Re_axpkxpz Im_axpkxpz];
            Re_axpkyz=real(axpkyz_mat(end)); Im_axpkyz=imag(axpkyz_mat(end)); axpkyz_mat_mat=[Re_axpkyz Im_axpkyz];
            Re_axpkypz=real(axpkypz_mat(end)); Im_axpkypz=imag(axpkypz_mat(end)); axpkypz_mat_mat=[Re_axpkypz Im_axpkypz];
            
            Re_aykd=real(aykd_mat(end)); Im_aykd=imag(aykd_mat(end)); aykd_mat_mat=[Re_aykd Im_aykd];
            Re_aykp=real(aykp_mat(end)); Im_aykp=imag(aykp_mat(end)); aykp_mat_mat=[Re_aykp Im_aykp];
            Re_aykxz=real(aykxz_mat(end)); Im_aykxz=imag(aykxz_mat(end)); aykxz_mat_mat=[Re_aykxz Im_aykxz];
            Re_aykxpz=real(aykxpz_mat(end)); Im_aykxpz=imag(aykxpz_mat(end)); aykxpz_mat_mat=[Re_aykxpz Im_aykxpz];
            Re_aykyz=real(aykyz_mat(end)); Im_aykyz=imag(aykyz_mat(end)); aykyz_mat_mat=[Re_aykyz Im_aykyz];
            Re_aykypz=real(aykypz_mat(end)); Im_aykypz=imag(aykypz_mat(end)); aykypz_mat_mat=[Re_aykypz Im_aykypz];
            
            Re_aypkd=real(aypkd_mat(end)); Im_aypkd=imag(aypkd_mat(end)); aypkd_mat_mat=[Re_aypkd Im_aypkd];
            Re_aypkp=real(aypkp_mat(end)); Im_aypkp=imag(aypkp_mat(end)); aypkp_mat_mat=[Re_aypkp Im_aypkp];
            Re_aypkxz=real(aypkxz_mat(end)); Im_aypkxz=imag(aypkxz_mat(end)); aypkxz_mat_mat=[Re_aypkxz Im_aypkxz];
            Re_aypkxpz=real(aypkxpz_mat(end)); Im_aypkxpz=imag(aypkxpz_mat(end)); aypkxpz_mat_mat=[Re_aypkxpz Im_aypkxpz];
            Re_aypkyz=real(aypkyz_mat(end)); Im_aypkyz=imag(aypkyz_mat(end)); aypkyz_mat_mat=[Re_aypkyz Im_aypkyz];
            Re_aypkypz=real(aypkypz_mat(end)); Im_aypkypz=imag(aypkypz_mat(end)); aypkypz_mat_mat=[Re_aypkypz Im_aypkypz];
        
            if (iscan==1)
                file_gkd=fopen(sprintf('gkd_mat_%02d.bin',file_index_write),'a');
                file_gkp=fopen(sprintf('gkp_mat_%02d.bin',file_index_write),'a');
                file_gkxz=fopen(sprintf('gkxz_mat_%02d.bin',file_index_write),'a');
                file_gkxpz=fopen(sprintf('gkxpz_mat_%02d.bin',file_index_write),'a');
                file_gkyz=fopen(sprintf('gkyz_mat_%02d.bin',file_index_write),'a');
                file_gkypz=fopen(sprintf('gkypz_mat_%02d.bin',file_index_write),'a');
                file_ekd=fopen(sprintf('ekd_mat_%02d.bin',file_index_write),'a');
                file_ekp=fopen(sprintf('ekp_mat_%02d.bin',file_index_write),'a');
                file_ekxz=fopen(sprintf('ekxz_mat_%02d.bin',file_index_write),'a');
                file_ekxpz=fopen(sprintf('ekxpz_mat_%02d.bin',file_index_write),'a');
                file_ekyz=fopen(sprintf('ekyz_mat_%02d.bin',file_index_write),'a');
                file_ekypz=fopen(sprintf('ekypz_mat_%02d.bin',file_index_write),'a');
                file_axkd=fopen(sprintf('axkd_mat_%02d.bin',file_index_write),'a');
                file_axkp=fopen(sprintf('axkp_mat_%02d.bin',file_index_write),'a');
                file_axkxz=fopen(sprintf('axkxz_mat_%02d.bin',file_index_write),'a');
                file_axkxpz=fopen(sprintf('axkxpz_mat_%02d.bin',file_index_write),'a');
                file_axkyz=fopen(sprintf('axkyz_mat_%02d.bin',file_index_write),'a');
                file_axkypz=fopen(sprintf('axkypz_mat_%02d.bin',file_index_write),'a');
                file_axpkd=fopen(sprintf('axpkd_mat_%02d.bin',file_index_write),'a');
                file_axpkp=fopen(sprintf('axpkp_mat_%02d.bin',file_index_write),'a');
                file_axpkxz=fopen(sprintf('axpkxz_mat_%02d.bin',file_index_write),'a');
                file_axpkxpz=fopen(sprintf('axpkxpz_mat_%02d.bin',file_index_write),'a');
                file_axpkyz=fopen(sprintf('axpkyz_mat_%02d.bin',file_index_write),'a');
                file_axpkypz=fopen(sprintf('axpkypz_mat_%02d.bin',file_index_write),'a');
                file_aykd=fopen(sprintf('aykd_mat_%02d.bin',file_index_write),'a');
                file_aykp=fopen(sprintf('aykp_mat_%02d.bin',file_index_write),'a');
                file_aykxz=fopen(sprintf('aykxz_mat_%02d.bin',file_index_write),'a');
                file_aykxpz=fopen(sprintf('aykxpz_mat_%02d.bin',file_index_write),'a');
                file_aykyz=fopen(sprintf('aykyz_mat_%02d.bin',file_index_write),'a');
                file_aykypz=fopen(sprintf('aykypz_mat_%02d.bin',file_index_write),'a');
                file_aypkd=fopen(sprintf('aypkd_mat_%02d.bin',file_index_write),'a');
                file_aypkp=fopen(sprintf('aypkp_mat_%02d.bin',file_index_write),'a');
                file_aypkxz=fopen(sprintf('aypkxz_mat_%02d.bin',file_index_write),'a');
                file_aypkxpz=fopen(sprintf('aypkxpz_mat_%02d.bin',file_index_write),'a');
                file_aypkyz=fopen(sprintf('aypkyz_mat_%02d.bin',file_index_write),'a');
                file_aypkypz=fopen(sprintf('aypkypz_mat_%02d.bin',file_index_write),'a');
            end
            
            fwrite(file_gkd,gkd_mat_mat,'double');
            fwrite(file_gkp,gkp_mat_mat,'double');
            fwrite(file_gkxz,gkxz_mat_mat,'double');
            fwrite(file_gkxpz,gkxpz_mat_mat,'double');
            fwrite(file_gkyz,gkyz_mat_mat,'double');
            fwrite(file_gkypz,gkypz_mat_mat,'double');
            
            fwrite(file_ekd,ekd_mat_mat,'double');
            fwrite(file_ekp,ekp_mat_mat,'double');
            fwrite(file_ekxz,ekxz_mat_mat,'double');
            fwrite(file_ekxpz,ekxpz_mat_mat,'double');
            fwrite(file_ekyz,ekyz_mat_mat,'double');
            fwrite(file_ekypz,ekypz_mat_mat,'double');
            
            fwrite(file_axkd,axkd_mat_mat,'double');
            fwrite(file_axkp,axkp_mat_mat,'double');
            fwrite(file_axkxz,axkxz_mat_mat,'double');
            fwrite(file_axkxpz,axkxpz_mat_mat,'double');
            fwrite(file_axkyz,axkyz_mat_mat,'double');
            fwrite(file_axkypz,axkypz_mat_mat,'double');
            
            fwrite(file_axpkd,axpkd_mat_mat,'double');
            fwrite(file_axpkp,axpkp_mat_mat,'double');
            fwrite(file_axpkxz,axpkxz_mat_mat,'double');
            fwrite(file_axpkxpz,axpkxpz_mat_mat,'double');
            fwrite(file_axpkyz,axpkyz_mat_mat,'double');
            fwrite(file_axpkypz,axpkypz_mat_mat,'double');
            
            fwrite(file_aykd,aykd_mat_mat,'double');
            fwrite(file_aykp,aykp_mat_mat,'double');
            fwrite(file_aykxz,aykxz_mat_mat,'double');
            fwrite(file_aykxpz,aykxpz_mat_mat,'double');
            fwrite(file_aykyz,aykyz_mat_mat,'double');
            fwrite(file_aykypz,aykypz_mat_mat,'double');
            
            fwrite(file_aypkd,aypkd_mat_mat,'double');
            fwrite(file_aypkp,aypkp_mat_mat,'double');
            fwrite(file_aypkxz,aypkxz_mat_mat,'double');
            fwrite(file_aypkxpz,aypkxpz_mat_mat,'double');
            fwrite(file_aypkyz,aypkyz_mat_mat,'double');
            fwrite(file_aypkypz,aypkypz_mat_mat,'double');
        
            if (iscan==scan_num)
                fclose('all');
            end
        end
    end   
    %}
end

%
% calculate final slice energy spread (SES), need to integrate over gain
% spectrum
if (scan_num>=20 && iEnergy_gain_calc==1 && iSES==1)
    if (iIBS==0); s_IBS=s; end
    %SES_fin=SES_calc(lambda_array,ekd_mat_fin_abs,Gf_11);     % unit: MeV
    SES_fin=SES_calc(lambda_array,ekd_mat_fin_abs,Gf_11,Gc_total);
end

if (checkbox31==1)
    generate_CSR_LSC_param;
end
%}
%--------------------- Post-processing -----------------------------------%
abs_method=1; %take abs, then sum
%abs_method=2; %sum, then take abs
if (scan_num==1)
    if (iplot_gain_func==1)
        if (iEnergy_gain_calc==0)
            figure(201); set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_11,'k-','linewidth',5);      xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        elseif (iEnergy_gain_calc==1)
            %figure(202); [hAx,hLine1,hLine2]=plotyy(s/100,G_11,s/100,G_22); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='-'; hLine2.LineStyle='-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{11}(s), density-to-density'); ylabel(hAx(2),'G_{22}(s), energy-to-density');
            %figure(203); [hAx,hLine1,hLine2]=plotyy(s/100,G_31,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='-'; hLine2.LineStyle='-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{31}(s), density-to-energy'); ylabel(hAx(2),'G_{44}(s), energy-to-energy');
            %figure(202); [hAx,hLine1,hLine2]=plotyy(s/100,G_11,s/100,G_22); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40); hLine1.LineStyle='-'; hLine2.LineStyle='-'; xlabel('s (m)'); ylabel(hAx(1),'G_{11}(s), density-to-density'); ylabel(hAx(2),'G_{22}(s), energy-to-density');
            %figure(203); [hAx,hLine1,hLine2]=plotyy(s/100,G_31,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40); hLine1.LineStyle='-'; hLine2.LineStyle='-'; xlabel('s (m)'); ylabel(hAx(1),'G_{31}(s), density-to-energy'); ylabel(hAx(2),'G_{44}(s), energy-to-energy');
            %
            %figure(204); [hAx,hLine1,hLine2]=plotyy(s/100,G_42,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='-'; hLine2.LineStyle='-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{42}(s)'); ylabel(hAx(2),'G_{44}(s)');
            figure(201); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b^d(s) (a.u.)'); hold on; grid off; axis('tight');
            %figure(206); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(G0_k),'b-','linewidth',5);   xlabel('s (m)'); ylabel('b_0^d(s) (a.u.)'); hold on; grid off; axis('tight');
            %figure(207); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b^e(s) (a.u.)'); hold on; grid off; axis('tight');
            %figure(208); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(G0_kp),'b-','linewidth',5);  xlabel('s (m)'); ylabel('b_0^e(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(209); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('p^d(s) (a.u.)'); hold on; grid off; axis('tight');
            %figure(210); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(E0_k),'b-','linewidth',5);   xlabel('s (m)'); ylabel('p_0^d(s) (a.u.)'); hold on; grid off; axis('tight');
            %figure(211); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('p^e(s) (a.u.)'); hold on; grid off; axis('tight');
            %figure(212); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(E0_kp),'b-','linewidth',5);  xlabel('s (m)'); ylabel('p_0^e(s) (a.u.)'); hold on; grid off; axis('tight');
            if (concatenate==1)
                if (iEnergy_gain_calc==1 && iTransverse_gain_calc==0)
                    if (abs_method==1)
                        figure(401); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat)+abs(gkp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b_k(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(402); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat)+abs(ekp_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('p_k(s) (a.u.)'); hold on; grid off; axis('tight');
                    end
                    if (abs_method==2)
                        figure(401); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat+gkp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b_k(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(402); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat+ekp_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('p_k(s) (a.u.)'); hold on; grid off; axis('tight');
                    end
                elseif (iEnergy_gain_calc==1 && iTransverse_gain_calc==1)
                    if (abs_method==1)
                        figure(401); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat)+abs(gkp_mat)+abs(gkxz_mat)+abs(gkxpz_mat)+abs(gkyz_mat)+abs(gkypz_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b_k(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(402); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat)+abs(ekp_mat)+abs(ekxz_mat)+abs(ekxpz_mat)+abs(ekyz_mat)+abs(ekypz_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('p_k(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(403); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(axkd_mat)+abs(axkp_mat)+abs(axkxz_mat)+abs(axkxpz_mat)+abs(axkyz_mat)+abs(axkypz_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(xz)}(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(404); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(axpkd_mat)+abs(axpkp_mat)+abs(axpkxz_mat)+abs(axpkxpz_mat)+abs(axpkyz_mat)+abs(axpkypz_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(x\primez)}(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(405); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(aykd_mat)+abs(aykp_mat)+abs(aykxz_mat)+abs(aykxpz_mat)+abs(aykyz_mat)+abs(aykypz_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(yz)}(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(406); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(aypkd_mat)+abs(aypkp_mat)+abs(aypkxz_mat)+abs(aypkxpz_mat)+abs(aypkyz_mat)+abs(aypkypz_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(y\primez)}(s) (a.u.)'); hold on; grid off; axis('tight');                        
                    end
                    if (abs_method==2)
                        figure(401); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat+gkp_mat+gkxz_mat+gkxpz_mat+gkyz_mat+gkypz_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b_k(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(402); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat+ekp_mat+ekxz_mat+ekxpz_mat+ekyz_mat+ekypz_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('p_k(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(403); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(axkd_mat+axkp_mat+axkxz_mat+axkxpz_mat+axkyz_mat+axkypz_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(xz)}(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(404); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(axpkd_mat+axpkp_mat+axpkxz_mat+axpkxpz_mat+axpkyz_mat+axpkypz_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(x\primez)}(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(405); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(aykd_mat+aykp_mat+aykxz_mat+aykxpz_mat+aykyz_mat+aykypz_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(yz)}(s) (a.u.)'); hold on; grid off; axis('tight');
                        figure(406); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(aypkd_mat+aypkp_mat+aypkxz_mat+aypkxpz_mat+aypkyz_mat+aypkypz_mat),'b-','linewidth',5);xlabel('s (m)'); ylabel('a_k^{(y\primez)}(s) (a.u.)'); hold on; grid off; axis('tight');                        
                    end
                end
            end
        end
    end
    if (Derbenev==1)
        figure(213); set(gca,'FontSize',40,'linewidth',5); plot(s(1:end-1)/100,Derbenev_ratio(1,:),'k','linewidth',5); xlabel('s (m)'); ylabel('Derbenev ratio \kappa=\sigma_x/\lambda^{2/3}|\rho|^{1/3}'); hold on; grid off; axis('tight');
    end
    if (first_harmonic_notification==1 && iEnergy_gain_calc==1)
        figure(214); set(gca,'FontSize',40,'linewidth',5); plot(s/100,LHS_01,'r','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
        figure(214); set(gca,'FontSize',40,'linewidth',5); plot(s/100,LHS_02,'g','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
        figure(214); set(gca,'FontSize',40,'linewidth',5); plot(s/100,LHS_03,'b','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
        figure(214); set(gca,'FontSize',40,'linewidth',5); plot(s/100,RHS,'k','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
    end
end

if (scan_num>=2)
    if (iplot_gain_spec==1)
        if (iEnergy_gain_calc==0)
            figure(300); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
            figure(314); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,max_G_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('max. gain'); hold on; grid off; axis('tight');
        elseif (iEnergy_gain_calc==1)
            %
            figure(9);  set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,11}, density-to-density'); hold on; grid off; axis('tight');
            figure(314); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,max_G_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('max. gain'); hold on; grid off; axis('tight');
            %figure(10); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_22(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,22}, energy-to-density'); hold on; grid off; axis('tight');
            figure(11); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_31(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,31},density-to-energy'); hold on; grid off; axis('tight');
            %figure(12); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_33(1:scan_num),'ko-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,33}'); hold on; grid off; axis('tight');
            %figure(13); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_42(1:scan_num),'ko-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,42}'); hold on; grid off; axis('tight');
            %figure(14); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_44(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,44}, energy-to-energy'); hold on; grid off; axis('tight');            
            %}
            %{
            figure(301); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,gkd_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|g_k^{d}|, density-to-density'); hold on; grid off; axis('tight');            
            figure(302); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,gkp_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|g_k^{e}|, energy-to-density'); hold on; grid off; axis('tight');            
            figure(303); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,ekd_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|p_k^{d}|, density-to-energy'); hold on; grid off; axis('tight');            
            figure(304); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,ekp_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|p_k^{e}|, energy-to-energy'); hold on; grid off; axis('tight');            
            %}
        end
        
        if (surfplot_gain_spec_func==1)
            figure(305); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,abs(G_total)); view(-38,32); xlabel('\lambda (\mum)','FontSize',40); ylabel('s (m)','FontSize',40); zlabel('G^d(s,\lambda)','FontSize',40); axis('tight'); 
            shading('interp'); h=colorbar; set(h,'FontSize',30);
            if (iEnergy_gain_calc==1)
                figure(306); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,abs(E_total)); view(-38,32); xlabel('\lambda (\mum)','FontSize',40); ylabel('s (m)','FontSize',40); zlabel('E^d(s,\lambda)','FontSize',40); axis('tight'); 
                shading('interp'); h=colorbar; set(h,'FontSize',30);
            end
        end
        
    end
    if (Derbenev==1)
        figure(307); set(gca,'FontSize',40,'linewidth',5); surf(s(1:end-1)/100,lambda_array*10^4,Derbenev_ratio); view(-140,50); xlabel('s (m)','FontSize',40); ylabel('\lambda (\mum)','FontSize',40); zlabel('Derbenev ratio \kappa=\sigma_x/\lambda^{2/3}|\rho|^{1/3}','FontSize',40); axis('tight'); shading interp;
    end
    if (first_harmonic_notification==1 && iEnergy_gain_calc==1)
        figure(308); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,LHS_01); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), d-e'); hold on; grid off; axis('tight'); shading interp;
        figure(309); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,LHS_02); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), e-e'); hold on; grid off; axis('tight'); shading interp;
        figure(310); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,LHS_03); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), total'); hold on; grid off; axis('tight'); shading interp;
        figure(311); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,RHS); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), d-e'); hold on; grid off; axis('tight'); shading interp;
        figure(312); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,RHS); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), e-e'); hold on; grid off; axis('tight'); shading interp;
        figure(313); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,RHS); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), total'); hold on; grid off; axis('tight'); shading interp;
    end
end

%if (isave==1) save('workplace.mat'); end
save('workplace.mat');
toc;

if (save_SDDS==1); gen_SDDS_MBI_data; end

fprintf('=============================================================================================== \n');
fprintf('Thanks for using volterra_mat.  Please cite the following reference in your publications: \n');
fprintf('C.-Y. Tsai, D. Douglas, R. Li, and C. Tennant, "Linear microbunching analysis for recirculation \n');
fprintf('machines" Phys. Rev. Accel. Beams 19, 114401 (2016). \n');
fprintf('If you use a modified version, please indicate this in all publications. \n');
fprintf('=============================================================================================== \n');

%volterra_plotter;
%elegant_plotter;

