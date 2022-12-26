
format long

global energy find_TWLA RF_info;

eegamma=energy*1000/0.511;

iplot=0;
prompt='Traveling-wave linac (TWLA) or standing-wave linac (RFCA) element in the lattice? (1:TWLA, 2:RFCA)';
find_TWLA=input(prompt);

%----------------- load "RF_elements_info.o" -----------------------------%
RF_info_gen_script;
% format (s,L,VOLT,PHASE,FREQ)
filename2='RF_elements_info.o';
delimiterIn=' '; headerlinesIn=0;
RF_info=importdata(filename2,delimiterIn,headerlinesIn);
[r,c]=size(RF_info);
if (mod(r,2)==0)
    NumberOfRFElements=r/2;
else
    fprintf('warning: RF information may be insufficient...\n');
end
%-------------------------------------------------------------------------%

%{
%----------------- load "pcentral_function.o" ----------------------------%
%pcentral_func_correction_script;
% format (s,pCentral=beta*gamma)
filename3='pcentral_function.o';
delimiterIn=' '; headerlinesIn=0;
pcentral_info=importdata(filename3,delimiterIn,headerlinesIn);
pcentral_s=pcentral_info(:,1);
pcentral_betagamma=pcentral_info(:,2);
%-------------------------------------------------------------------------%

if (find_TWLA==1) %TWLA
%******** additional contribution of R56 due to drfit section ********% 
drift_correction=1;
if (drift_correction==1)
    egamma_vec=interp1(pcentral_s,pcentral_betagamma,s_ele);
	R56_drift=s_ele./egamma_vec.^3;
	transport(:,31)=transport(:,31)+R56_drift;
end
egamma_vec_old=egamma_vec;
%*********************************************************************%


for m=1:1:NumberOfRFElements
    
    TWLA_s_start=RF_info(2*m-1,1);
    TWLA_length=RF_info(2*m-1,2);
    TWLA_s_end=TWLA_s_start+TWLA_length;
    s_mesh_num=TWLA_length/0.01;
    TWLA_freq=RF_info(2*m-1,5);
    TWLA_psi0_deg=RF_info(2*m-1,4);
    TWLA_E0=RF_info(2*m-1,3);
    
    egamma0=interp1(pcentral_s,pcentral_betagamma,TWLA_s_start);
    
    M_ini(1,1)=interp1(s_ele,transport(:,2),TWLA_s_start);
    M_ini(1,2)=interp1(s_ele,transport(:,3),TWLA_s_start);
    M_ini(1,3)=interp1(s_ele,transport(:,4),TWLA_s_start);
    M_ini(1,4)=interp1(s_ele,transport(:,5),TWLA_s_start);
    M_ini(1,5)=interp1(s_ele,transport(:,6),TWLA_s_start);
    M_ini(1,6)=interp1(s_ele,transport(:,7),TWLA_s_start);
    M_ini(2,1)=interp1(s_ele,transport(:,8),TWLA_s_start);
    M_ini(2,2)=interp1(s_ele,transport(:,9),TWLA_s_start);
    M_ini(2,3)=interp1(s_ele,transport(:,10),TWLA_s_start);
    M_ini(2,4)=interp1(s_ele,transport(:,11),TWLA_s_start);
    M_ini(2,5)=interp1(s_ele,transport(:,12),TWLA_s_start);
    M_ini(2,6)=interp1(s_ele,transport(:,13),TWLA_s_start);
    M_ini(3,1)=interp1(s_ele,transport(:,14),TWLA_s_start);
    M_ini(3,2)=interp1(s_ele,transport(:,15),TWLA_s_start);
    M_ini(3,3)=interp1(s_ele,transport(:,16),TWLA_s_start);
    M_ini(3,4)=interp1(s_ele,transport(:,17),TWLA_s_start);
    M_ini(3,5)=interp1(s_ele,transport(:,18),TWLA_s_start);
    M_ini(3,6)=interp1(s_ele,transport(:,19),TWLA_s_start);
    M_ini(4,1)=interp1(s_ele,transport(:,20),TWLA_s_start);
    M_ini(4,2)=interp1(s_ele,transport(:,21),TWLA_s_start);
    M_ini(4,3)=interp1(s_ele,transport(:,22),TWLA_s_start);
    M_ini(4,4)=interp1(s_ele,transport(:,23),TWLA_s_start);
    M_ini(4,5)=interp1(s_ele,transport(:,24),TWLA_s_start);
    M_ini(4,6)=interp1(s_ele,transport(:,25),TWLA_s_start);
    M_ini(5,1)=interp1(s_ele,transport(:,26),TWLA_s_start);
    M_ini(5,2)=interp1(s_ele,transport(:,27),TWLA_s_start);
    M_ini(5,3)=interp1(s_ele,transport(:,28),TWLA_s_start);
    M_ini(5,4)=interp1(s_ele,transport(:,29),TWLA_s_start);
    M_ini(5,5)=interp1(s_ele,transport(:,30),TWLA_s_start);
    M_ini(5,6)=interp1(s_ele,transport(:,31),TWLA_s_start);
    M_ini(6,1)=interp1(s_ele,transport(:,32),TWLA_s_start);
    M_ini(6,2)=interp1(s_ele,transport(:,33),TWLA_s_start);
    M_ini(6,3)=interp1(s_ele,transport(:,34),TWLA_s_start);
    M_ini(6,4)=interp1(s_ele,transport(:,35),TWLA_s_start);
    M_ini(6,5)=interp1(s_ele,transport(:,36),TWLA_s_start);
    M_ini(6,6)=interp1(s_ele,transport(:,37),TWLA_s_start);
    
    total_Y=RF_gap_mat_6x6_linear(egamma0,M_ini,TWLA_s_start,TWLA_length,TWLA_freq,TWLA_psi0_deg,TWLA_E0,s_mesh_num,m);
    
    M_fin=total_Y(end,1:37);
    
    %---------------------------------------------------------------------%
    % search the row index (of lattice_transport_functions.o) of which the starting position of TWLA begins
    [~,ind_start]=min(abs(s_ele-TWLA_s_start));
    [~,ind_end]=min(abs(s_ele-TWLA_s_end));
    
    % replace the original transport functions by the above 6x6 RF gap transfer map at "the" specific index
    transport_tmp(1:ind_start-1,:)=transport(1:ind_start-1,:);
    transport_tmp(ind_start:ind_start+s_mesh_num-1,:)=total_Y(:,1:37);
    %egamma_vec_tmp(1:ind_start-1,:)=egamma_vec(1:ind_start-1,:);
    %egamma_vec_tmp(ind_start:ind_start+s_mesh_num-1,:)=total_Y(:,38);
    
    % update 6x6 transfer map for all downstream elements
    [r,~]=size(transport);
    transport_diff=M_fin-transport(ind_end,:);
    transport_tmp(ind_start+s_mesh_num:ind_start+s_mesh_num+r-ind_end,:)=transport(ind_end:end,:)+repmat(transport_diff,r-ind_end+1,1);
    transport=transport_tmp;
    
    
    % remove duplicate s to avoid error in interpolation
    [r,~]=size(transport);
    %tmp(1,:)=transport(1,:); n=2; q=0;
    tmp=transport(1,:); n=2; q=0;
    for p=2:1:r
        if (transport(p-1,1)==transport(p,1))
            q=q+1;
        else
            tmp(n,:)=transport(p,:); 
            n=n+1;
        end
    end
    clear transport;
    transport=tmp;
    s_ele=transport(:,1);
    %---------------------------------------------------------------------%
    
end

if (iplot==1)
figure(1001);
subplot(2,2,1); plot(transport(:,1),transport(:,30),'b'); hold on; xlabel('s (m)'); ylabel('R_{55} (m)');
subplot(2,2,2); plot(transport(:,1),transport(:,31),'b'); hold on; xlabel('s (m)'); ylabel('R_{56}');
subplot(2,2,3); plot(transport(:,1),transport(:,36),'b'); hold on; xlabel('s (m)'); ylabel('R_{65} (1/m)');
subplot(2,2,4); plot(transport(:,1),transport(:,37),'b'); hold on; xlabel('s (m)'); ylabel('R_{66}');
end

% base conversion from (x,px,y,py,z,Delta_gamma) to (x,px,y,py,z,delta)
egamma_vec=interp1(pcentral_s,pcentral_betagamma,s_ele);
transport(:,7)=transport(:,7).*egamma_vec;
transport(:,13)=transport(:,13).*egamma_vec;
transport(:,19)=transport(:,19).*egamma_vec;
transport(:,25)=transport(:,25).*egamma_vec;
transport(:,31)=transport(:,31).*egamma_vec;
transport(:,32)=transport(:,32)./egamma_vec;
transport(:,33)=transport(:,33)./egamma_vec;
transport(:,34)=transport(:,34)./egamma_vec;
transport(:,35)=transport(:,35)./egamma_vec;
transport(:,36)=transport(:,36)./egamma_vec;


% minus-sign conversion back to ELEGANT format, now z<0 for bunch head
transport(:,6)=-transport(:,6);
transport(:,12)=-transport(:,12);
transport(:,18)=-transport(:,18);
transport(:,24)=-transport(:,24);
transport(:,36)=-transport(:,36);
transport(:,26)=-transport(:,26);
transport(:,27)=-transport(:,27);
transport(:,28)=-transport(:,28);
transport(:,29)=-transport(:,29);
transport(:,31)=-transport(:,31);

fileID1='lattice_transport_functions.o';
dlmwrite(fileID1,transport,'delimiter',' ','precision','%10.10f');

pcentral_data=[s_ele egamma_vec];
fileID2='pcentral_function.o';
dlmwrite(fileID2,pcentral_data,'delimiter',' ','precision','%10.10f');

if (iplot==1)
figure(1002);
subplot(2,2,1); plot(transport(:,1),transport(:,30),'b'); hold on; xlabel('s (m)'); ylabel('R_{55} (m)');
subplot(2,2,2); plot(transport(:,1),transport(:,31),'b'); hold on; xlabel('s (m)'); ylabel('R_{56}');
subplot(2,2,3); plot(transport(:,1),transport(:,36),'b'); hold on; xlabel('s (m)'); ylabel('R_{65} (1/m)');
subplot(2,2,4); plot(transport(:,1),transport(:,37),'b'); hold on; xlabel('s (m)'); ylabel('R_{66}');

figure(1003); plot(transport(:,1),egamma_vec,'b'); hold on; xlabel('s (m)'); ylabel('\gamma_r');

figure(1004);
subplot(2,2,1); plot(transport(:,1),transport(:,26),'b'); hold on; xlabel('s (m)'); ylabel('R_{51}');
subplot(2,2,2); plot(transport(:,1),transport(:,27),'b'); hold on; xlabel('s (m)'); ylabel('R_{52} (m)');
subplot(2,2,3); plot(transport(:,1),transport(:,28),'b'); hold on; xlabel('s (m)'); ylabel('R_{53}');
subplot(2,2,4); plot(transport(:,1),transport(:,29),'b'); hold on; xlabel('s (m)'); ylabel('R_{54} (m)');

figure(1005);
subplot(2,2,1); plot(transport(:,1),transport(:,2),'black'); hold on; xlabel('s (m)'); ylabel('R_{11}');
subplot(2,2,2); plot(transport(:,1),transport(:,3),'black'); hold on; xlabel('s (m)'); ylabel('R_{12} (m)');
subplot(2,2,3); plot(transport(:,1),transport(:,8),'black'); hold on; xlabel('s (m)'); ylabel('R_{21}');
subplot(2,2,4); plot(transport(:,1),transport(:,9),'black'); hold on; xlabel('s (m)'); ylabel('R_{22} (m)');
end
end
%}


