function [s_m alphax_m alphay_m betax_m betay_m psix_m psiy_m etax_m etay_m]=dipole_info_Twiss_new(flag)

%clear all; clc;
format long;

filename01='elements_geom.o';    % format (ElementType,s)
filename02='Twiss_lattice.o';    % format (s,alphax,alphay,betax,betay,psix,psiy)
delimiterIn01='\t'; delimiterIn02=' '; headerlinesIn=0;
dipole_location=importdata(filename01,delimiterIn01,headerlinesIn);
dipole_Twiss=importdata(filename02,delimiterIn02,headerlinesIn);

[r,~]=size(dipole_location.data);


%============ select dipole elements and read dipole information =========%
n=0; p=0; q=0; w=0;
for m=1:1:r
    
    if ((strcmp(dipole_location.textdata(m),'CSRCSBEND')==1)||...
        (strcmp(dipole_location.textdata(m),'CSBEND')==1)||...
        (strcmp(dipole_location.textdata(m),'KSBEND')==1)||...
        (strcmp(dipole_location.textdata(m),'NISBEND')==1)||...
        (strcmp(dipole_location.textdata(m),'RBEN')==1)||...
        (strcmp(dipole_location.textdata(m),'SBEN')==1)||...
        (strcmp(dipole_location.textdata(m),'TUBEND')==1))
        w=w+1;
        if (w==1)
            n=n+1;
            bend_ent_row_index(n)=m-1;
        end
    else
        if (n~=0 && w~=0)
            p=p+1;
            bend_ext_row_index(p)=m-1;
            w=0;
        end
    end
    
end
%=========================================================================%
if (n~=0)
[rr,cc]=size(bend_ent_row_index); n=cc;

bend_ent_ext_row_index(1:2:2*cc-1)=bend_ent_row_index;
bend_ent_ext_row_index(2:2:2*cc)=bend_ext_row_index;

%-------------------------------------------------------------------------%
[rrr,ccc]=size(bend_ent_ext_row_index);
if (ccc ~= 2*cc)
    fprintf('warning: something wrong in lattice_dipole_info.o...\n');
    fprintf('remember to add at least one watch point or marker at the end of lattice...\n');
    fprintf('press any key to continue or Ctrl+C to terminate...\n');
    pause;
end
%-------------------------------------------------------------------------%

dipole_ent_s=dipole_location.data(bend_ent_row_index,1);
dipole_exit_s=dipole_location.data(bend_ext_row_index,1);

%============ format conversion ==========================================%
dipole_s(1:2:2*n-1)=dipole_ent_s;
dipole_s(2:2:2*n)=dipole_exit_s;

shortest=min(dipole_exit_s-dipole_ent_s)*100;
%=========================================================================%

%============ summary and output corrected lattice file ==================%
fprintf('generating the file lattice_dipole_info.o...'); 
fprintf('done... \n');
fprintf('total number of %d dipoles detected... \n',cc);
fprintf('shortest dipole length in the lattice = %f cm \n',shortest);

else
    dipole_mat=[0;1e50;1e50]';
    fprintf('no dipole detected...\n'); shortest=1.0; % cm
end

%============== calculate middle point of each dipole ====================%
s_m=0.5*(dipole_ent_s+dipole_exit_s);
%=========================================================================%

%================= extract dipole Twiss parameters =======================%
dipole_s_ref=dipole_Twiss(:,1);

%============ check if duplicate s coordinates are present ===============%
[r,~]=size(dipole_s_ref);
dipole_s_ref_new(1,:)=dipole_s_ref(1,:); n=2; q=0;
for m=2:1:r
    if (dipole_s_ref(m-1,1)==dipole_s_ref(m,1))
        q=q+1;
    else
        dipole_s_ref_new(n,:)=dipole_s_ref(m,:);
        dipole_Twiss_new(n,:)=dipole_Twiss(m,:);
        n=n+1;
    end
end
dipole_s_ref=dipole_s_ref_new;
%=========================================================================%

dipole_alphax=dipole_Twiss_new(:,2);
dipole_alphay=dipole_Twiss_new(:,3);
dipole_betax= dipole_Twiss_new(:,4);
dipole_betay= dipole_Twiss_new(:,5);
dipole_psix=  dipole_Twiss_new(:,6);
dipole_psiy=  dipole_Twiss_new(:,7);
dipole_etax=  dipole_Twiss_new(:,8);
dipole_etay=  dipole_Twiss_new(:,9);

alphax_m=interp1(dipole_s_ref_new,dipole_alphax,s_m);
alphay_m=interp1(dipole_s_ref_new,dipole_alphay,s_m);
betax_m= interp1(dipole_s_ref_new,dipole_betax,s_m);
betay_m= interp1(dipole_s_ref_new,dipole_betay,s_m);
psix_m=  interp1(dipole_s_ref_new,dipole_psix,s_m);
psiy_m=  interp1(dipole_s_ref_new,dipole_psiy,s_m);
etax_m=  interp1(dipole_s_ref_new,dipole_etax,s_m);
etay_m=  interp1(dipole_s_ref_new,dipole_etay,s_m);

alphax_i=interp1(dipole_s_ref_new,dipole_alphax,dipole_ent_s);
alphay_i=interp1(dipole_s_ref_new,dipole_alphay,dipole_ent_s);
betax_i= interp1(dipole_s_ref_new,dipole_betax,dipole_ent_s);
betay_i= interp1(dipole_s_ref_new,dipole_betay,dipole_ent_s);
psix_i=  interp1(dipole_s_ref_new,dipole_psix,dipole_ent_s);
psiy_i=  interp1(dipole_s_ref_new,dipole_psiy,dipole_ent_s);
etax_i=  interp1(dipole_s_ref_new,dipole_etax,dipole_ent_s);
etay_i=  interp1(dipole_s_ref_new,dipole_etay,dipole_ent_s);

alphax_f=interp1(dipole_s_ref_new,dipole_alphax,dipole_exit_s);
alphay_f=interp1(dipole_s_ref_new,dipole_alphay,dipole_exit_s);
betax_f= interp1(dipole_s_ref_new,dipole_betax,dipole_exit_s);
betay_f= interp1(dipole_s_ref_new,dipole_betay,dipole_exit_s);
psix_f=  interp1(dipole_s_ref_new,dipole_psix,dipole_exit_s);
psiy_f=  interp1(dipole_s_ref_new,dipole_psiy,dipole_exit_s);
etax_f=  interp1(dipole_s_ref_new,dipole_etax,dipole_exit_s);
etay_f=  interp1(dipole_s_ref_new,dipole_etay,dipole_exit_s);
%=========================================================================%
%{
figure(96); set(gca,'FontSize',40,'linewidth',5); plot(1:1:cc,alphax_m,'ro--'); xlabel('index of dipoles'); ylabel('\alpha function'); hold on;
figure(97); set(gca,'FontSize',40,'linewidth',5); plot(1:1:cc,betax_m,'bo--');  xlabel('index of dipoles'); ylabel('\beta function (m)'); hold on;
figure(98); set(gca,'FontSize',40,'linewidth',5); plot(1.5:1:(cc-0.5),diff(psix_m/pi),'blacko--'); xlabel('index of dipoles'); ylabel('N-N phase difference ({\pi})'); hold on;
 
figure(99);
[ax,h1,h2]=plotyy(1.5:1:(cc-0.5),diff(psix_m/pi),1:1:cc,etax_m); xlabel('index of dipoles'); ylabel('N-N phase difference ({\pi})'); hold on;
set(ax(1),'Position',[0.13 0.11+0.05 0.775 0.815]);
set(ax(2),'Position',[0.13 0.11+0.05 0.775 0.815]);
%}
%
figure(96); set(gca,'FontSize',40,'linewidth',5); errorbar(1:1:cc,0.5*abs(alphax_i-alphax_f),'ro--'); xlabel('index of dipoles'); ylabel('\alpha function'); hold on;
figure(97); set(gca,'FontSize',40,'linewidth',5); errorbar(1:1:cc,0.5*abs(betax_i-betax_f),'bo--');  xlabel('index of dipoles'); ylabel('\beta function (m)'); hold on;
figure(98); set(gca,'FontSize',40,'linewidth',5); errorbar(1.5:1:(cc-0.5),diff(psix_m/pi),'blacko--'); xlabel('index of dipoles'); ylabel('N-N phase difference ({\pi})'); hold on;
 
figure(99);
[ax,h1,h2]=plotyy(1.5:1:(cc-0.5),diff(psix_m/pi),1:1:cc,etax_m); xlabel('index of dipoles'); ylabel('N-N phase difference ({\pi})'); hold on;
set(ax(1),'Position',[0.13 0.11+0.05 0.775 0.815]);
set(ax(2),'Position',[0.13 0.11+0.05 0.775 0.815]);
%}
%{
%================ plot s_m_i vs. s_m_j vs. diff(si-sj) ===================%
[r,~]=size(psix_m);
horizontal_index=1:1:r;
vertical_index=1:1:r;

for m=2:1:r
    for n=1:1:r
        sisj_map(m,n)=psix_m(m)-psix_m(n);
    end
end

figure(4); surf(mod((sisj_map/pi),1)); view(0,90); colorbar;
%=========================================================================%
%}

end
