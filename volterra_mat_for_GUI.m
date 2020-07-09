%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                         %
%% This program calculates microbunching gain based on linearized Vlasov   %
%% equation (in mathematical form of Volterra integral equation), including%
%% collective effects such as coherent synchrotron radiation (CSR),        %
%% longitudinal space charge (LSC) and linac geometric impedances.         %
%% The impedance models are in analytical forms. CSR is assumed 1-D, with  %
%% ultrarelativistic (for high-energy bean) and non-ultrarelativistic      %
%% (for low-energy beam) versions. There are two options for steady-state  %
%% CSR models: free-space model and parallel-plate shielding model. The    %
%% transient CSR models, including entrance transient and exit transient   %
%% (or, CSR drift), are based on Saldin's treatment. The LSC models have   %
%% various versions, for different transverse beam distributions and       %
%% various cross sections of surrounding pipes. See Yingjie Li's thesis    %
%% (Michigan State. U., 2015) or step-by-step guide for more details.      %
%%                                                                         %
%% As of v3.3 and later, we have included both density and energy          %
%% microbunching. For more details, see C.-Y. Tsai's Ph.D. thesis          %
%% (Virginia Tech, 2017) or a published paper (NIMA 940, 462-474, 2019).   %
%%                                                                         %
%% The input requires physics parameters and numerical parameters. The     %
%% former are read from ELEGANT, so users need to prepare two ELEGANT files%
%% (*.ele and *.lte) for a certain transport line, and specify a set of    %
%% output formats (*.twi,*.mag,*.sig,*.cen,*.flr,*.mat,*.param, etc).      %
%% The latter are provided by users, through GUI.                          %
%%                                                                         %
%% This program calculates the gain function G(s), if a specific modulation%
%% wavelength is given or final gain spectrum Gf, if a series of modulation%
%% wavelengths are scanned.                                                %
%%                                                                         %
%% For users interested in the original theoretical formulation, we refer  %
%% them to the paper by Heifets, Stupakov, and Krinsky, PRST-AB 5, 064401  %
%% (2002), and Huang and Kim, PRST-AB 5, 074401 (2002) and their errata.   %
%%                                                                         %
%% Benchmarking of this program against ELEGANT was documented in          %
%% JLAB-TN-14-016. There was a step-by-step guide to using this program,   %
%% which was documented in JLAB-TN-15-019. Benchmarking against energy     %
%% modulation induced microbunching can be found in JLAB-TN-16-022.        %
%% Although this code has been benchmarked against many lattice examples,  %
%% it may still have bugs. Users should be careful about unusual outputs.  %
%% The step-by-step guide is recently updated (July, 2020).                %
%%                                                                         %
%% Note: The options with finite-bunch or Gaussian bunch distribution, the %
%%       transverse-longitudinal microbunching, iterative calculation, and %
%%       validity check of first-harmonics are not available in this       %
%%       version. The option with inclusion of RF element is not included  %
%%       in this version and is available in v4.2.1_acc_r version.         %
%%                                                                         %
%% Program written by Cheng-Ying Tsai. For any comments/questions, please  % 
%% send to jcytsai@hust.edu.cn                                             %
%%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
%------------------ Global variable declaration --------------------------%
global egamma s_ele R11_ele R12_ele R13_ele R14_ele R21_ele R22_ele R23_ele R24_ele R31_ele R32_ele R33_ele R34_ele R41_ele R42_ele R43_ele R44_ele R16_ele R26_ele R36_ele R46_ele R51_ele R52_ele R53_ele R54_ele R55_ele R56_ele C_ele;
global rhox rhoy k_wave emit_norm_x emit_norm_y emitx emity alphax0 alphay0 betax0 betay0 gammax0 gammay0 sigma_delta n_1k0 e_1k0 re nb I_b chirp start_pos end_pos I_A const_A mesh_num;
global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac issCSRpp;
global sig_x_ele sig_y_ele egamma_vec find_TWLA full_pipe_height round_pipe_radius;
global LSC_model ssCSR_model first;
global Nb sigma_z0;

format long
const_A=3^(-1/3)*gamma(2/3)*(1i*sqrt(3)-1);
c_speed=2.99792458e10;      % unit: cm/sec
charge=1.6e-19;             % unit: Coulomb
first=1; % signal indicating shielding CSR model is used (see kernel_mod.m)

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
HH = findobj(gcf,'Tag','iplot_lattice');   iplot_lattice = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iplot_gain_func'); iplot_gain_func = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iplot_gain_spec'); iplot_gain_spec = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iplot_energy_mod');iplot_energy_mod = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','isave');           isave = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iiterative');      iiterative = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','ilast');           ilast = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iCSR_ss');         iCSR_ss = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iCSR_tr');         iCSR_tr = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iCSR_drift');      iCSR_drift = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iLSC');            iLSC = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','ilinac');          ilinac = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iGaussian');       iGaussian = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','iEnergy_mod_calc');iEnergy_mod_calc = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iEnergy_gain_calc');iEnergy_gain_calc = str2num(get(HH,'String'));
%HH = findobj(gcf,'Tag','idm_analysis');    idm_analysis = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','surfplot_gain_spec_func');surfplot_gain_spec_func = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','issCSRpp');        issCSRpp = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','ssCSR_model');     ssCSR_model = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','LSC_model');       LSC_model = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','round_pipe_radius');round_pipe_radius = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iquilt_plot');     iquilt_plot = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iplot_I_b');       iplot_I_b = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iRF_ele');         iRF_ele = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','full_pipe_height');full_pipe_height = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','Derbenev');        Derbenev = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iTransverse_gain_calc');iTransverse_gain_calc = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','first_harmonic_notification');first_harmonic_notification = str2num(get(HH,'String'));

%--------------------- Initialize beam parameters ------------------------%
egamma=energy*1000/0.511;                           % Lorentz factor
I_A=17045;                                          % Alfven current in Amp
re=2.81794*10^(-13);                                % classical radius in cm
nb=I_b/(re*I_A);                                    % line density (cm^-1)
Nb=2.35*sigma_z0*I_b/(c_speed*charge);              % total number of physical particles
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
iRF_ele=0;
if (iRF_ele==0); find_TWLA=0; end                   % this option is disabled in this version

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
%enx_s_ele=beam_sigma_mat(:,6);                     % enx_s in m
%eny_s_ele=beam_sigma_mat(:,7);                     % eny_s in m

%C_ele=sig_s_ele(1)./sig_s_ele;

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

%--------------------- Start of the main program -------------------------%
for iscan=1:1:scan_num
    
    n_1k0=1;
    e_1k0=1;
    
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
    
%----------------------Solve Volterra integral equation-------------------%
    
    s=linspace(start_pos,end_pos,mesh_num);
        
    G0_k=g0k_mat(s);                                % density bunching due to density modulation
    if (iEnergy_gain_calc==1)
        G0_kp=g0kp_mat(s);                          % density bunching due to energy  modulation
        E0_k=p0k_mat(s);                            % energy  bunching due to density modulation
        E0_kp=p0kp_mat(s);                          % energy  bunching due to energy  modulation
    end

    q=1;                                            % dipole index
    K_mat=zeros(mesh_num);
    M_minus_L_mat=zeros(mesh_num);
    
    if (Derbenev==1)
        for p=1:1:(mesh_num-1)
            Derbenev_ratio(iscan,p)=interp1(s_ele,sig_x_ele,s(p))/(lambda^(2/3)*abs(auxr(s(p)))^(1/3));
        end
    end
    
    if ((iCSR_drift==1) || (iLSC==1) || (ilinac==1))
    for p=1:1:(mesh_num-1)
    	s_q=s(p+1:mesh_num);
        if (p==1)
            K_mat(p+1:mesh_num,p)=0.5*kernel_mod(s_q,s(p));  
            if (iEnergy_gain_calc==1)
                M_minus_L_mat(p+1:mesh_num,p)=0.5*(kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p)));
            end
        else
            K_mat(p+1:mesh_num,p)=kernel_mod(s_q,s(p));
            if (iEnergy_gain_calc==1)
                M_minus_L_mat(p+1:mesh_num,p)=kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p));
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
                end
            else
                K_mat(p+1:mesh_num,p)=kernel_mod(s_q,s(p));
                if (iEnergy_gain_calc==1)
                    M_minus_L_mat(p+1:mesh_num,p)=kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p));
                end
            end
        end
        
    end
    end
        
    % -------------- self-consistent (direct) solution ------------------ %
    
    if (iEnergy_gain_calc==0)
        gkd_mat=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_k);
        G_c_11=gkd_mat/G0_k(1);
        G_11=abs(G_c_11);
        Gf_11(iscan)=G_11(end);
        G_total(:,iscan)=G_11;
    elseif (iEnergy_gain_calc==1 && iTransverse_gain_calc==0)
        tot_col_vec_new(1:mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_k);
        tot_col_vec_new(1*mesh_num+1:2*mesh_num)=(eye(mesh_num)-mesh_size*K_mat)\transpose(G0_kp);
        tot_col_vec_new(2*mesh_num+1:3*mesh_num)=mesh_size*M_minus_L_mat*((eye(mesh_num)-mesh_size*K_mat)\transpose(G0_k))+transpose(E0_k);
        tot_col_vec_new(3*mesh_num+1:4*mesh_num)=mesh_size*M_minus_L_mat*((eye(mesh_num)-mesh_size*K_mat)\transpose(G0_kp))+transpose(E0_kp);
    end
    
    if (iEnergy_gain_calc==1)       
        gkd_mat=tot_col_vec_new(1:mesh_num);
        gkp_mat=tot_col_vec_new(mesh_num+1:2*mesh_num); 
        ekd_mat=tot_col_vec_new(2*mesh_num+1:3*mesh_num);
        ekp_mat=tot_col_vec_new(3*mesh_num+1:4*mesh_num);
        
        G_c_11=gkd_mat/G0_k(1); G_11=abs(G_c_11); Gf_11(iscan)=G_11(end);
        %G_c_22=gkp_mat/G0_kp(1);G_22=abs(G_c_22); Gf_22(iscan)=G_22(end);
        G_c_22=gkp_mat/e_1k0;   G_22=abs(G_c_22); Gf_22(iscan)=G_22(end);
        G_c_31=ekd_mat/G0_k(1); G_31=abs(G_c_31); Gf_31(iscan)=G_31(end);
        %G_c_33=ekd_mat/E0_k(1); G_33=abs(G_c_33); Gf_33(iscan)=G_33(end);
        %G_c_42=ekp_mat/G0_kp(1);G_42=abs(G_c_42); Gf_42(iscan)=G_42(end);
        %G_c_44=ekp_mat/E0_kp(1);G_44=abs(G_c_44); Gf_44(iscan)=G_44(end);
        G_c_44=ekp_mat/e_1k0;   G_44=abs(G_c_44); Gf_44(iscan)=G_44(end);
        
        G_total(:,iscan)=abs(G_c_11+G_c_22);
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
    end
    
    if (iEnergy_gain_calc==1)
    	fprintf('iscan = %4d, lambda = %6.2f um, d-d gain = %5.2f, d-e gain = %5.2f, e-d gain = %5.2f, e-e gain = %5.2f \n',iscan,lambda*10^4,Gf_11(iscan),Gf_31(iscan),Gf_22(iscan),Gf_44(iscan));
    else
    	fprintf('iscan = %4d, lambda = %6.2f um, density gain = %5.2f \n',iscan,lambda*10^4,Gf_11(iscan));
    end
end


%--------------------- Post-processing -----------------------------------%
abs_method=1; %take abs, then sum
%abs_method=2; %sum, then take abs
if (scan_num==1)
    if (iplot_gain_func==1)
        if (iEnergy_gain_calc==0)
            figure(201); set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_11,'k-','linewidth',5);      xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        elseif (iEnergy_gain_calc==1)
            %figure(202); [hAx,hLine1,hLine2]=plotyy(s/100,G_11,s/100,G_22); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{11}(s), density-to-density'); ylabel(hAx(2),'G_{22}(s), energy-to-density');
            %figure(203); [hAx,hLine1,hLine2]=plotyy(s/100,G_31,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{31}(s), density-to-energy'); ylabel(hAx(2),'G_{44}(s), energy-to-energy');
            %figure(204); [hAx,hLine1,hLine2]=plotyy(s/100,G_42,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{42}(s)'); ylabel(hAx(2),'G_{44}(s)');
            figure(205); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b^d(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(206); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(G0_k),'b-','linewidth',5);   xlabel('s (m)'); ylabel('b_0^d(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(207); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b^e(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(208); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(G0_kp),'b-','linewidth',5);  xlabel('s (m)'); ylabel('b_0^e(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(209); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('p^d(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(210); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(E0_k),'b-','linewidth',5);   xlabel('s (m)'); ylabel('p_0^d(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(211); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('p^e(s) (a.u.)'); hold on; grid off; axis('tight');
            figure(212); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(E0_kp),'b-','linewidth',5);  xlabel('s (m)'); ylabel('p_0^e(s) (a.u.)'); hold on; grid off; axis('tight');
        end
    end
    if (Derbenev==1)
        figure(213); set(gca,'FontSize',40,'linewidth',5); plot(s(1:end-1)/100,Derbenev_ratio(1,:),'k','linewidth',5); xlabel('s (m)'); ylabel('Derbenev ratio \kappa=\sigma_x/\lambda^{2/3}|\rho|^{1/3}'); hold on; grid off; axis('tight');
    end
end

if (scan_num>=2)
    if (iplot_gain_spec==1)
        if (iEnergy_gain_calc==0)
            figure(300); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
        elseif (iEnergy_gain_calc==1)
            %{
            figure(9);  set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,11}, density-to-density'); hold on; grid off; axis('tight');
            figure(10); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_22(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,22}, energy-to-density'); hold on; grid off; axis('tight');
            figure(11); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_31(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,31},density-to-energy'); hold on; grid off; axis('tight');
            %figure(12); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_33(1:scan_num),'ko-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,33}'); hold on; grid off; axis('tight');
            %figure(13); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_42(1:scan_num),'ko-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,42}'); hold on; grid off; axis('tight');
            figure(14); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_44(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,44}, energy-to-energy'); hold on; grid off; axis('tight');            
            %}
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
end
%if (isave==1) save('workplace.mat'); end
save('workplace.mat');
toc;

fprintf('=============================================================================================== \n');
fprintf('Thanks for using volterra_mat. Please cite the following reference in your publications: \n');
fprintf('C.-Y. Tsai, D. Douglas, R. Li, and C. Tennant, "Linear microbunching analysis for \n');
fprintf('recirculation machines," Phys. Rev. Accel. Beams 19, 114401 (2016). \n');
fprintf('=============================================================================================== \n');

%volterra_plotter;
%elegant_plotter;

