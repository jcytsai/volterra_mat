%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This program calculates microbunching gain based on linearized Vlasov   %
% equation (or, Volterra integral equation in mathematical form), with    %
% collective effects such as coherent synchrotron radiation (CSR),        %
% longitudinal space charge(LSC) and linac geometric impedances considered%
% The impedance models are in analytical form. CSR is assumed 1-D, with   %
% ultra-relativistic (usually used for high-energy case) and relativistic %
% (can be used for low-energy case) versions. There are two options for   %
% CSR models; the free space model and parallel-plate shielding model.    %
% The LSC models have various versions, for different transverse beam     %
% distributions and various shapes (cross sections) of surrounding pipes. %
%                                                                         %
% The required input information includes physical parameters and         %
% numerical parameters. The former is read from ELEGANT, so users need to %
% prepare a set of ELEGANT input files (*.ele and *.lte) for a certain    %
% lattice, and specify a set of ELEGANT output files (*.twi, *.mag, *.sig,%
% *.cen, *.flr, *.mat, *.param, etc) The latter must be provided by users,%
% through the interactive message or GUI.                                 %
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
% pleae refer to our IPAC'16 work.                                        %
%                                                                         %
% Benchmarking of this program against ELEGANT was documented in          %
% JLAB-TN-14-016. There is a step-by-step guide to using this program,    %
% which was documented in JLAB-TN-15-017.                                 %
%                                                                         %
% Details about impedance models used in the code can be found in:        %
% C. -Y. Tsai et al., IPAC'15 (MODMA028)                                  %
%                                                                         %
% Note: since v3.3, we have removed iterative (or staged) gain calculation%
% from the program.                                                       %
%                                                                         %
% Program written by Cheng-Ying Tsai, jcytsai@vt.edu                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
%------------------ Global variable declaration --------------------------%
global egamma s_ele R16_ele R36_ele R51_ele R52_ele R53_ele R54_ele R55_ele R56_ele C_ele rhox rhoy k_wave emit_norm_x emit_norm_y emitx emity alphax0 alphay0 betax0 betay0 sigma_delta n_1k0 e_1k0 re nb I_b chirp start_pos end_pos I_A const_A mesh_num;
%global dipole_s C_factor lambda_start01 lambda_end01 scan_num01 mesh_num;
global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac issCSRpp;
global sig_x_ele sig_y_ele egamma_vec find_TWLA full_pipe_height;
global LSC_model ssCSR_model first;
global dipole_control_seq;

format long
const_A=3^(-1/3)*gamma(2/3)*(1i*sqrt(3)-1);
first=1; % signal indicating shielding CSR model is used (see kernel_mod.m)
Derbenev=1;
first_harmonic_notification=1;

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
HH = findobj(gcf,'Tag','isave');           isave = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iiterative');      iiterative = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','ilast');           ilast = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iCSR_ss');         iCSR_ss = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iCSR_tr');         iCSR_tr = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iCSR_drift');      iCSR_drift = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iLSC');            iLSC = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','ilinac');          ilinac = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iGaussian');       iGaussian = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iEnergy_mod_calc');iEnergy_mod_calc = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iEnergy_gain_calc');iEnergy_gain_calc = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','idm_analysis');    idm_analysis = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','surfplot_gain_spec_func');surfplot_gain_spec_func = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','issCSRpp');        issCSRpp = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','ssCSR_model');     ssCSR_model = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','LSC_model');       LSC_model = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iquilt_plot');     iquilt_plot = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iplot_I_b');       iplot_I_b = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','iRF_ele');         iRF_ele = str2num(get(HH,'String'));
HH = findobj(gcf,'Tag','full_pipe_height');full_pipe_height = str2num(get(HH,'String'));


%--------------------- Initialize beam parameters ------------------------%
egamma=energy*1000/0.511;                           % Lorentz factor
I_A=17045;                                          % Alfven current in Amp
re=2.81794*10^(-13);                                % classical radius in cm
nb=I_b/(re*I_A);                                    % line density (cm^-1)
emit_norm_x=emit_norm_x*10^(-4);                    % normalized horizontal emittance in cm
emitx=emit_norm_x/egamma;                           % geometric  horizontal emittance in cm
emit_norm_y=emit_norm_y*10^(-4);                    % normalized vertical   emittance in cm
emity=emit_norm_y/egamma;                           % geometric  vertical   emittance in cm
betax0=betax0*10^2;                                 % initial betax in cm
betay0=betay0*10^2;                                 % initial betay in cm
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
    fprintf('ignore the message if care of dipole action control is not made...\n');
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
R16_ele=100*transport(:,7);                         % dispersion Dx(s)
R36_ele=100*transport(:,13);                        % dispersion Dy(s)
R51_ele=transport(:,26);                            % R51 unitless
R52_ele=100*transport(:,27);                        % R52 in cm
R53_ele=transport(:,28);                            % R53 unitless
R54_ele=100*transport(:,29);                        % R54 in cm
R55_ele=transport(:,30);                            % R55 unitless
R56_ele=-100*transport(:,31);                       % R56 in cm

%********* conversion of R51: from ELEGANT to HSK format *****************%
R51_ele=-R51_ele+R52_ele*(alphax0/betax0);
R52_ele=-R52_ele;
if (betay0 > 0)
    R53_ele=-R53_ele+R54_ele*(alphay0/betay0);
    R54_ele=-R54_ele;
end

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

%************************** Plot lattice transport functions *************%
if (iplot_lattice==1)
	figure(1); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R16_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{16} (cm)');      grid off; hold on; axis('tight');
	figure(2); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R36_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{36} (cm)');      grid off; hold on; axis('tight');
	figure(3); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R51_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{51}');           grid off; hold on; axis('tight');
	figure(4); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R52_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{52} (cm)');      grid off; hold on; axis('tight');
	figure(5); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R53_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{53}');           grid off; hold on; axis('tight');
	figure(6); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R54_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{54} (cm)');      grid off; hold on; axis('tight');
    figure(7); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,R56_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('R_{56} (cm)');      grid off; hold on; axis('tight');
    figure(8); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele,'b-','linewidth',5); xlabel('s (m)'); ylabel('Compression factor'); grid off; hold on; axis('tight');
    %figure(8); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele*I_b,'b-','linewidth',5); xlabel('s (m)'); ylabel('Peak current (A)'); grid off; hold on; axis('tight');
end
%*************************************************************************%

%************************** Plot beam current evolution ******************%
if (iplot_I_b==1); figure(888); set(gca,'FontSize',40,'linewidth',5); plot(s_ele/100,C_ele*I_b,'b-','linewidth',5); xlabel('s (m)'); ylabel('Peak current (A)'); grid off; hold on; axis('tight'); end
%*************************************************************************%

%******* Read transverse rms beam size function for LSC calculation ******%
filename='beam_sigma_mat_functions.o';
% foramt (s,sigma_x,sigma_y,sigma_s,sigma_delta,emit_norm_x,emit_norm_y)
beam_sigma_mat_tmp=importdata(filename,delimiterIn,headerlinesIn);    
beam_sigma_mat=interp1(beam_sigma_mat_tmp(:,1),beam_sigma_mat_tmp,s_ele/100);
sig_x_ele=100*beam_sigma_mat(:,2);                  % sigma_x in cm
sig_y_ele=100*beam_sigma_mat(:,3);                  % sigma_y in cm
%sig_s_ele=100*beam_sigma_mat(:,4);                 % sigma_s in cm
%sdelta_s_ele= beam_sigma_mat(:,5);
%enx_s_ele=beam_sigma_mat(:,6);                     % enx_s in m
%eny_s_ele=beam_sigma_mat(:,7);                     % eny_s in m
%*************************************************************************%

%************************** Plot quilt pattern R56(s'->s) ****************%
if (iquilt_plot==1); quilt_plot; end
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

%--------------------- Option for concatenation study --------------------%
concatenate=0;
if (concatenate==1)
    [tmp01,tmp02,tmp03,tmp04]=read_tot_col_vec_4x4(2,scan_num);
    n_1k0_array=tmp01+tmp02;
    e_1k0_array=tmp03+tmp04;
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
    
    if (concatenate==0)
        n_1k0=1.0;                                  % a constant
        e_1k0=1.0;                                  % a constant
    else
        n_1k0=n_1k0_array(iscan);
        e_1k0=e_1k0_array(iscan);
    end
%----------------------Solve Volterra integral equation-------------------%
    
    s=linspace(start_pos,end_pos,mesh_num);
        
    G0_k=g0k_mat(s);                                % density bunching due to density modulation without CSR
    if (iEnergy_gain_calc==1)
        G0_kp=g0kp_mat(s);                          % density bunching due to energy  modulation without CSR
        E0_k=p0k_mat(s);                            % energy  bunching due to density modulation without CSR
        E0_kp=p0kp_mat(s);                          % energy  bunching due to energy  modulation without CSR
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
                %M_minus_L_mat(p+1:mesh_num,p)=0.5*(-kernel_mod_L(s_q,s(p)));
            end
        else
            K_mat(p+1:mesh_num,p)=kernel_mod(s_q,s(p));
            if (iEnergy_gain_calc==1)
                M_minus_L_mat(p+1:mesh_num,p)=kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p));
                %M_minus_L_mat(p+1:mesh_num,p)=-kernel_mod_L(s_q,s(p));
            end
        end
    end
    
    else
    
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
                    %M_minus_L_mat(p+1:mesh_num,p)=0.5*(-kernel_mod_L(s_q,s(p)));
                end
            else
                K_mat(p+1:mesh_num,p)=kernel_mod(s_q,s(p));
                if (iEnergy_gain_calc==1)
                    M_minus_L_mat(p+1:mesh_num,p)=kernel_mod_M(s_q,s(p))-kernel_mod_L(s_q,s(p));
                    %M_minus_L_mat(p+1:mesh_num,p)=-kernel_mod_L(s_q,s(p));
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
    elseif (iEnergy_gain_calc==1)
        one_minus_K_mat_inv=inv(eye(mesh_num)-mesh_size*K_mat);
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
        %{
        % read tot_col_vec_old from external file
        prompt='read tot_col_vec_old from external files? (1:yes,0:no)';
        ans=input(prompt);
        if (ans==1)
            tot_col_vec_old=read_tot_col_vec(mesh_num,2);
        end
        %}
        tot_col_vec_new=MBI_kernel_mat*tot_col_vec_old;
    end
    
    if (iEnergy_gain_calc==1)
        gkd_mat=tot_col_vec_new(1:mesh_num);
        gkp_mat=tot_col_vec_new(mesh_num+1:2*mesh_num); 
        ekd_mat=tot_col_vec_new(2*mesh_num+1:3*mesh_num);
        ekp_mat=tot_col_vec_new(3*mesh_num+1:end);
        
        G_c_11=gkd_mat/G0_k(1); G_11=abs(G_c_11); Gf_11(iscan)=G_11(end);
        %G_c_22=gkp_mat/G0_kp(1);G_22=abs(G_c_22); Gf_22(iscan)=G_22(end);
        G_c_22=gkp_mat/e_1k0;        G_22=abs(G_c_22); Gf_22(iscan)=G_22(end);
        G_c_31=ekd_mat/G0_k(1); G_31=abs(G_c_31); Gf_31(iscan)=G_31(end);
        %G_c_33=ekd_mat/E0_k(1); G_33=abs(G_c_33); Gf_33(iscan)=G_33(end);
        %G_c_42=ekp_mat/G0_kp(1);G_42=abs(G_c_42); Gf_42(iscan)=G_42(end);
        %G_c_44=ekp_mat/E0_kp(1);G_44=abs(G_c_44); Gf_44(iscan)=G_44(end);
        G_c_44=ekp_mat/e_1k0;        G_44=abs(G_c_44); Gf_44(iscan)=G_44(end);
        
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
        
        if (first_harmonic_notification==1)
            %if (iscan==1)
                LHS_01(:,iscan)=abs(interp1(s_ele,R56_ele,s))'.*abs(ekd_mat);          % unit: cm
                LHS_02(:,iscan)=abs(interp1(s_ele,R56_ele,s))'.*abs(ekp_mat);          % unit: cm
                LHS_03(:,iscan)=abs(interp1(s_ele,R56_ele,s))'.*abs(ekd_mat+ekp_mat);  % unit: cm
                RHS(:,iscan)=lambda./interp1(s_ele,C_ele,s);                           % unit: cm
            %end
        end
        
    end
    
    if (iEnergy_gain_calc==1)
    	fprintf('iscan = %d, lambda = %f um, emit_norm_x = %f um, emit_norm_y = %f um, dp/p = %f, d-d gain = %f, d-e gain = %f, e-d gain = %f, e-e gain = %f \n',iscan,lambda*10^4,emit_norm_x*10^4,emit_norm_y*10^4,sigma_delta,Gf_11(iscan),Gf_31(iscan),Gf_22(iscan),Gf_44(iscan));
    else
    	fprintf('iscan = %d, lambda = %f um, emit_norm_x = %f um, emit_norm_y = %f um, dp/p = %f, density gain = %f \n',iscan,lambda*10^4,emit_norm_x*10^4,emit_norm_y*10^4,sigma_delta,Gf_11(iscan));
    end
    
    if (scan_num>1 && iEnergy_gain_calc==1)   
        ans=0;
        if (ans==1)
            Re_gkd=real(gkd_mat(end)); Im_gkd=imag(gkd_mat(end)); gkd_mat_mat=[Re_gkd Im_gkd];
            Re_gkp=real(gkp_mat(end)); Im_gkp=imag(gkp_mat(end)); gkp_mat_mat=[Re_gkp Im_gkp];
            Re_ekd=real(ekd_mat(end)); Im_ekd=imag(ekd_mat(end)); ekd_mat_mat=[Re_ekd Im_ekd];
            Re_ekp=real(ekp_mat(end)); Im_ekp=imag(ekp_mat(end)); ekp_mat_mat=[Re_ekp Im_ekp];
        
            if (iscan==1)
            file_gkd=fopen('gkd_mat_02.bin','a');
            file_gkp=fopen('gkp_mat_02.bin','a');
            file_ekd=fopen('ekd_mat_02.bin','a');
            file_ekp=fopen('ekp_mat_02.bin','a');
            end
            
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
    end   
    
end


%--------------------- Post-processing -----------------------------------%

if (scan_num==1)
    if (iplot_gain_func==1)
        if (iEnergy_gain_calc==0)
            figure(88); set(gca,'FontSize',40,'linewidth',5); plot(s/100,G_11,'k-','linewidth',5);      xlabel('s (m)'); ylabel('G(s)'); hold on; grid off; axis('tight');
        elseif (iEnergy_gain_calc==1)
            figure(15); [hAx,hLine1,hLine2]=plotyy(s/100,G_11,s/100,G_22); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{11}(s), density-to-density'); ylabel(hAx(2),'G_{22}(s), energy-to-density');
            figure(16); [hAx,hLine1,hLine2]=plotyy(s/100,G_31,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{31}(s), density-to-energy'); ylabel(hAx(2),'G_{44}(s), energy-to-energy');
            %figure(17); [hAx,hLine1,hLine2]=plotyy(s/100,G_42,s/100,G_44); grid off; axis(hAx,'tight'); set(hAx,'FontSize',40,'linewidth',5); set(hLine1,'linewidth',5); set(hLine2,'linewidth',5); hLine1.LineStyle='r-'; hLine2.LineStyle='b-'; hLine1.Linewidth=5; hLine2.Linewidth=5; xlabel('s (m)'); ylabel(hAx(1),'G_{42}(s)'); ylabel(hAx(2),'G_{44}(s)');
            figure(18); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkd_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b^d(s)'); hold on; grid off; axis('tight');
            figure(18); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(G0_k),'b-','linewidth',5);   xlabel('s (m)'); ylabel('b^d(s)'); hold on; grid off; axis('tight');
            figure(19); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(gkp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('b^e(s)'); hold on; grid off; axis('tight');
            figure(19); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(G0_kp),'b-','linewidth',5);  xlabel('s (m)'); ylabel('b^e(s)'); hold on; grid off; axis('tight');
            figure(20); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekd_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('p^d(s)'); hold on; grid off; axis('tight');
            figure(20); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(E0_k),'b-','linewidth',5);   xlabel('s (m)'); ylabel('p^d(s)'); hold on; grid off; axis('tight');
            figure(22); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(ekp_mat),'r-','linewidth',5);xlabel('s (m)'); ylabel('p^e(s)'); hold on; grid off; axis('tight');
            figure(22); set(gca,'FontSize',40,'linewidth',5); plot(s/100,abs(E0_kp),'b-','linewidth',5);  xlabel('s (m)'); ylabel('p^e(s)'); hold on; grid off; axis('tight');
        end
    end
    if (Derbenev==1)
        figure(30); set(gca,'FontSize',40,'linewidth',5); plot(s(1:end-1)/100,Derbenev_ratio(1,:),'k','linewidth',5); xlabel('s (m)'); ylabel('Derbenev ratio \kappa=\sigma_x/\lambda^{2/3}|\rho|^{1/3}'); hold on; grid off; axis('tight');
    end
    if (first_harmonic_notification==1 && iEnergy_gain_calc==1)
        figure(31); set(gca,'FontSize',40,'linewidth',5); plot(s/100,LHS_01,'r','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
        figure(31); set(gca,'FontSize',40,'linewidth',5); plot(s/100,LHS_02,'g','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
        figure(31); set(gca,'FontSize',40,'linewidth',5); plot(s/100,LHS_03,'b','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
        figure(31); set(gca,'FontSize',40,'linewidth',5); plot(s/100,RHS,'k','linewidth',5); xlabel('s (m)'); ylabel('R_{56}(s)p(k;s), \lambda(s) (cm) (R:d-e,G:e-e,B:total)'); hold on; grid off; axis('tight');
    end
end

if (scan_num>=2)
    if (iplot_gain_spec==1)
        if (iEnergy_gain_calc==0)
            figure(9); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f}'); hold on; grid off; axis('tight');
        elseif (iEnergy_gain_calc==1)
            %{
            figure(9);  set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_11(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,11}, density-to-density'); hold on; grid off; axis('tight');
            figure(10); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_22(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,22}, energy-to-density'); hold on; grid off; axis('tight');
            figure(11); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_31(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,31},density-to-energy'); hold on; grid off; axis('tight');
            %figure(12); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_33(1:scan_num),'ko-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,33}'); hold on; grid off; axis('tight');
            %figure(13); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_42(1:scan_num),'ko-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,42}'); hold on; grid off; axis('tight');
            figure(14); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,Gf_44(1:scan_num),'b-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('G_{f,44}, energy-to-energy'); hold on; grid off; axis('tight');            
            %}
            figure(23); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,gkd_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|g_k^{d}|, density-to-density'); hold on; grid off; axis('tight');            
            figure(24); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,gkp_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|g_k^{e}|, energy-to-density'); hold on; grid off; axis('tight');            
            figure(25); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,ekd_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|p_k^{d}|, density-to-energy'); hold on; grid off; axis('tight');            
            figure(26); set(gca,'FontSize',40,'linewidth',5); plot(lambda_array(1:scan_num)*10^4,ekp_mat_fin_abs(1:scan_num),'r-','linewidth',5); xlabel('\lambda (\mum)'); ylabel('|p_k^{e}|, energy-to-energy'); hold on; grid off; axis('tight');            
            %}
        end
        
        if (surfplot_gain_spec_func==1)
            figure(27); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,abs(G_total)); view(-38,32); xlabel('\lambda (\mum)','FontSize',40); ylabel('s (m)','FontSize',40); zlabel('G^d(s,\lambda)','FontSize',40); axis('tight'); 
            shading('interp'); h=colorbar; set(h,'FontSize',30);
            figure(28); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,abs(E_total)); view(-38,32); xlabel('\lambda (\mum)','FontSize',40); ylabel('s (m)','FontSize',40); zlabel('E^d(s,\lambda)','FontSize',40); axis('tight'); 
            shading('interp'); h=colorbar; set(h,'FontSize',30);
        end
        
    end
    if (Derbenev==1)
        figure(30); set(gca,'FontSize',40,'linewidth',5); surf(s(1:end-1)/100,lambda_array*10^4,Derbenev_ratio); view(-140,50); xlabel('s (m)','FontSize',40); ylabel('\lambda (\mum)','FontSize',40); zlabel('Derbenev ratio \kappa=\sigma_x/\lambda^{2/3}|\rho|^{1/3}','FontSize',40); axis('tight'); shading interp;
    end
    if (first_harmonic_notification==1 && iEnergy_gain_calc==1)
        figure(31); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,LHS_01); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), d-e'); hold on; grid off; axis('tight'); shading interp;
        figure(32); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,LHS_02); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), e-e'); hold on; grid off; axis('tight'); shading interp;
        figure(33); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,LHS_03); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), total'); hold on; grid off; axis('tight'); shading interp;
        figure(31); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,RHS); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), d-e'); hold on; grid off; axis('tight'); shading interp;
        figure(32); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,RHS); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), e-e'); hold on; grid off; axis('tight'); shading interp;
        figure(33); set(gca,'FontSize',40,'linewidth',5); surf(lambda_array*10^4,s/100,RHS); ylabel('s (m)'); xlabel('\lambda (\mum)'); zlabel('R_{56}(s)p(k;s), \lambda(s) (cm), total'); hold on; grid off; axis('tight'); shading interp;
    end
end
%if (isave==1) save('workplace.mat'); end
save('workplace.mat');
toc;

fprintf('=============================================================================================== \n');
fprintf('Thanks for using volterra_mat.  Please cite the following reference in your publications: \n');
fprintf('C. -Y. Tsai, D. Douglas, R. Li, and C. Tennant, "Linear Vlasov Solver for Microbunching Gain \n');
fprintf('Estimation with Inclusion of CSR, LSC And Linac Geometric Impedances," FEL 2015 (MOP052) \n');
fprintf('If you use a modified version, please indicate this in all publications. \n');
fprintf('=============================================================================================== \n');

%volterra_plotter;

