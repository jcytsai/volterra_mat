%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subroutine is written to generate a linear map for a linac section,
% within which 6x6 transport functions are evaulated as a function of s.
% This subroutine would be called by trans_mat_trim.m and given the following
% input arguments:
% 
% Input: egamma0, initial beam energy
%        M_ini, initial 6x6 transfer matrix elements
%        TWLA_s_start, linac starting position
%        TWLA_length, length of linac section
%        TWLA_freq, RF frequency in the linac
%        TWLA_psi0_deg, initial RF phase w.r.t. peak voltage in deg
%        TWLA_E0, RF peak accelerating gradient or peak voltage (MeV/m)
%        s_mesh_num, number of s-mesh along the linac
%
% Output:(s,R11,R12,...,R65,R66.egamma), a matrix with size (s_mech_num)x38
%
% Program written by Cheng-Ying Tsai, jcytsai@vt.edu
% v1: April 6, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function func=RF_gap_mat_6x6_linear(egamma0,M_ini,TWLA_s_start,TWLA_length,TWLA_freq,TWLA_psi0_deg,TWLA_E0,s_mesh_num,index_flag)

format long

c_speed=2.99792458e8;
small=0.001;

iplot=0;
%--------- initial beam parameter and transport matrix setting -----------%
%egamma0=(250)/0.511;                % beam (kinetic) energy in MeV
%M_ini=eye(6);                             % 6x6 M(s0), should adapt from ELEGANT later
%-------------------------------------------------------------------------%

%--------- initial RF parameters setting----------------------------------%
%TWLA_s_start=0.0;                         % start position, should be given later an absolute s-coordinate
%TWLA_length=330;                          % unit: m
%TWLA_freq=2.855168171942223e+09;          % unit: Hz
%TWLA_psi0_deg=-43;                         % initial RF phase in deg
%TWLA_E0=20;                               % peak accelerating gradient in MeV/m
TWLA_end_pos=TWLA_s_start+TWLA_length;    % end position
%k_RF=2*pi*TWLA_freq/c_speed;
%-------------------------------------------------------------------------%

%s_mesh_num=400;
s=linspace(TWLA_s_start,TWLA_end_pos-small,s_mesh_num);

%--------- calculation of transport functions ----------------------------%

[egamma_vec,Y]=M_gen(TWLA_s_start,TWLA_length,egamma0,TWLA_freq,TWLA_psi0_deg,TWLA_E0,M_ini);

%-------------------------------------------------------------------------%

if (iplot==1)

figure(1);
subplot(2,2,1); plot(Y(:,1),Y(:,30),'b'); hold on; xlabel('s (m)'); ylabel('R_{55} (m)');
subplot(2,2,2); plot(Y(:,1),Y(:,31),'b'); hold on; xlabel('s (m)'); ylabel('R_{56}');
subplot(2,2,3); plot(Y(:,1),Y(:,36),'b'); hold on; xlabel('s (m)'); ylabel('R_{65} (m)');
subplot(2,2,4); plot(Y(:,1),Y(:,37),'b'); hold on; xlabel('s (m)'); ylabel('R_{66}');

figure(2);
plot(Y(:,1),egamma_vec,'r'); xlabel('s (m)'); ylabel('\gamma_r');

figure(3);
subplot(2,2,1); plot(Y(:,1),Y(:,26),'b'); hold on; xlabel('s (m)'); ylabel('R_{51}');
subplot(2,2,2); plot(Y(:,1),Y(:,27),'b'); hold on; xlabel('s (m)'); ylabel('R_{52} (m)');
subplot(2,2,3); plot(Y(:,1),Y(:,28),'b'); hold on; xlabel('s (m)'); ylabel('R_{53}');
subplot(2,2,4); plot(Y(:,1),Y(:,29),'b'); hold on; xlabel('s (m)'); ylabel('R_{54} (m)');

figure(4);
subplot(2,2,1); plot(Y(:,1),Y(:,2),'black'); hold on; xlabel('s (m)'); ylabel('R_{11}');
subplot(2,2,2); plot(Y(:,1),Y(:,3),'black'); hold on; xlabel('s (m)'); ylabel('R_{12} (m)');
subplot(2,2,3); plot(Y(:,1),Y(:,7),'black'); hold on; xlabel('s (m)'); ylabel('R_{21}');
subplot(2,2,4); plot(Y(:,1),Y(:,8),'black'); hold on; xlabel('s (m)'); ylabel('R_{22} (m)');
end

tmp_Y=interp1(Y(:,1),Y,s);
tmp_egamma_vec=interp1(Y(:,1),egamma_vec,s);

total_Y=[tmp_Y tmp_egamma_vec'];
%fileID1=sprintf('RF_transport_functions_%d.o',index_flag);
%dlmwrite(fileID1,total_Y,'delimiter',' ','precision','%.15f');

func=total_Y;





