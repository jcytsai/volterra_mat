% THIS PROGRAM SUGGESTS A SET OF NUMERICAL PARAMETERS FOR ELEGANT TRACKING
% BASED ON THE MICROBUNCHING GAIN SPECTRUM
% REFERENCE: M. BORLAND, OAG-TN-2005-027, C.-Y.TSAI & R. LI, JLAB-TN-14-016
% INPUT: Gf_11, lambda_array

%global Gf_max lambda_opt lambda_min num_bin_per_lambda num_ptc_per_bin sig_s_ini 
global num_kick num_bin num_ptc cutoff_01 cutoff_02 ini_dmod_amp ini_emod_amp

[val,ind]=max(abs(Gf_11));
Gf_max=val;
lambda_opt=lambda_array(ind);                     % unit: cm
lambda_min=lambda_array(1);                       % unit: cm
num_bin_per_lambda=10;
num_ptc_per_bin=3000;
sig_s_ini=sig_s_ele(1);

if (gen_CSR_LSC_param==1)                         % IDEALIZED CASE
    B_param=10*lambda_opt;                        % unit: cm
    M_param=lambda_min;                           % unit: cm
    num_bin=round(num_bin_per_lambda*B_param/M_param);
    num_ptc=round(num_ptc_per_bin*num_bin);
    cutoff_01=2/num_bin_per_lambda;
    cutoff_02=4/num_bin_per_lambda;
    num_kick=round(shortest/0.5);
    ini_dmod_amp=0.1/Gf_max;
    %ini_emod_amp=max(abs(Gf_31));
    
elseif (gen_CSR_LSC_param==2)                     % REALISTIC CASE
    B_param=sig_s_ini;                            % unit: cm
    M_param=lambda_min;
    
    if (B_param<4*M_param)
        fprintf('Bunch length comparable to minimum lambda...\n');
    end
    
    num_bin=round(num_bin_per_lambda*B_param/M_param);
    num_ptc=round(num_ptc_per_bin*num_bin);
    cutoff_01=2/num_bin_per_lambda;
    cutoff_02=4/num_bin_per_lambda;
    num_kick=round(shortest/0.5);
    ini_dmod_amp=0.1/Gf_max;
    %ini_emod_amp=max(abs(Gf_31));
end

GUI_CSR_LSC_param_gen;

%{
CSRCSBEND,L=0.2, ANGLE=0.0609, E2=0.0609, EDGE_ORDER=1,n_kicks=40,NONLINEAR=0,LINEARIZE=1,&
INTEGRATION_ORDER=4,BINS="nbins_csr",SG_HALFWIDTH=1,SGDERIV_HALFWIDTH=1,ISR="flag_isr",CSR="flag_csr",&
HIGH_FREQUENCY_CUTOFF0="f0_csr",HIGH_FREQUENCY_CUTOFF1="f1_csr",output_file="B_LH.01.CSR",OUTPUT_LAST_WAKE_ONLY=1

CSRDRIFT, DZ=0.01, CSR="flag_csr_drif", USE_STUPAKOV=1, LINEARIZE=1, LSC_INTERPOLATE="flag_lsc_drif",
LSC_BINS="nbins_lsc",LSC_LOW_FREQUENCY_CUTOFF0="f0_lsc",LSC_LOW_FREQUENCY_CUTOFF1="f1_lsc",
LSC_HIGH_FREQUENCY_CUTOFF0="f0_lsc",LSC_HIGH_FREQUENCY_CUTOFF1="f1_lsc", L=0.11

LSCDRIFT, L=1.9840, bins="nbins_lsc", smoothing=1, HIGH_FREQUENCY_CUTOFF0="f0_lsc",
HIGH_FREQUENCY_CUTOFF1="f1_lsc", lsc="flag_lsc"
%}
