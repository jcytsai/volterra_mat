function CSRSuppressionConditionCheck

% see also, dipole_info_Twiss.m

format long;

iplot=0;

% Read ELEGANT lattice transport functions, use MKS
load('lattice_transport_functions.o');
transport=lattice_transport_functions;
s_ele=transport(:,1);                           % Frenet-Serret s (in cm)
R16_ele=transport(:,7);                         % dispersion Dx(s)
R26_ele=transport(:,13);
R56_ele=-transport(:,31);                       % R56 in cm

% Read ELEGANT lattice dipole functions, use MKS
dipole_info_gen_script;

% Read ELEGANT lattice Twiss functions, use MKS
filename01='elements_geom.o';    % format (ElementType,s)
filename02='Twiss_lattice.o';    % format (s,alphax,alphay,betax,betay,psix,psiy)
delimiterIn01='\t'; delimiterIn02=' '; headerlinesIn=0;
dipole_location=importdata(filename01,delimiterIn01,headerlinesIn);
dipole_Twiss=importdata(filename02,delimiterIn02,headerlinesIn);

%============== calculate middle point of each dipole ====================%
s_m=0.5*(dipole_ent_s+dipole_exit_s);
%=========================================================================%

%================= extract dipole Twiss parameters =======================%
dipole_s_ref=dipole_Twiss(:,1);

%============ check if duplicate s coordinates are present ===============%
[r,~]=size(dipole_s_ref);
dipole_s_ref_new(1,:)=dipole_s_ref(1,:); n=1; q=0;
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

% i: initial, dipole_ent_s, alphax_i, ...
% m: middle, s_m, alphax_m, ...
% f: final, dipole_exit_s, alphax_f, ...

alphax_i=interp1(dipole_s_ref_new,dipole_alphax,dipole_ent_s,'linear','extrap');
alphay_i=interp1(dipole_s_ref_new,dipole_alphay,dipole_ent_s,'linear','extrap');
betax_i= interp1(dipole_s_ref_new,dipole_betax,dipole_ent_s,'linear','extrap');
betay_i= interp1(dipole_s_ref_new,dipole_betay,dipole_ent_s,'linear','extrap');
psix_i=  interp1(dipole_s_ref_new,dipole_psix,dipole_ent_s,'linear','extrap');
psiy_i=  interp1(dipole_s_ref_new,dipole_psiy,dipole_ent_s,'linear','extrap');
etax_i=  interp1(dipole_s_ref_new,dipole_etax,dipole_ent_s,'linear','extrap');
etay_i=  interp1(dipole_s_ref_new,dipole_etay,dipole_ent_s,'linear','extrap');

alphax_m=interp1(dipole_s_ref_new,dipole_alphax,s_m);
alphay_m=interp1(dipole_s_ref_new,dipole_alphay,s_m);
betax_m= interp1(dipole_s_ref_new,dipole_betax,s_m);
betay_m= interp1(dipole_s_ref_new,dipole_betay,s_m);
psix_m=  interp1(dipole_s_ref_new,dipole_psix,s_m);
psiy_m=  interp1(dipole_s_ref_new,dipole_psiy,s_m);
etax_m=  interp1(dipole_s_ref_new,dipole_etax,s_m);
etay_m=  interp1(dipole_s_ref_new,dipole_etay,s_m);

alphax_f=interp1(dipole_s_ref_new,dipole_alphax,dipole_exit_s);
alphay_f=interp1(dipole_s_ref_new,dipole_alphay,dipole_exit_s);
betax_f= interp1(dipole_s_ref_new,dipole_betax,dipole_exit_s);
betay_f= interp1(dipole_s_ref_new,dipole_betay,dipole_exit_s);
psix_f=  interp1(dipole_s_ref_new,dipole_psix,dipole_exit_s);
psiy_f=  interp1(dipole_s_ref_new,dipole_psiy,dipole_exit_s);
etax_f=  interp1(dipole_s_ref_new,dipole_etax,dipole_exit_s);
etay_f=  interp1(dipole_s_ref_new,dipole_etay,dipole_exit_s);

%{
%=========================================================================%
[r,~]=size(s_m);
scan_num=40;

flag01=0; flag02=0;
for m=1:1:(r-1)
    alphax_scan=linspace(0.1*alphax_i(m),10*alphax_i(m),scan_num);
    betax_scan=linspace(0.1*betax_i(m),5*betax_i(m),scan_num);
    Psi_scan=linspace(0,2*pi,20);
    s1=0.5*abs(dipole_ent_s(m)-dipole_exit_s(m));
    s2=0.5*abs(dipole_ent_s(m+1)-dipole_exit_s(m+1));
    Lb=2*s1;
    rho1=dipole_mat(2*m,2);
    rho2=dipole_mat(2*(m+1),2);
    
    if (rho1*rho2<0)
        s2=-s2;
        fprintf('%d-th to %d-th dipole reverse bend detected\n',m,m+1);
    end
    
    rho1=abs(rho1);
    rho2=abs(rho2);
    
    R16=abs(interp1(s_ele,R16_ele,dipole_exit_s(m+1))-interp1(s_ele,R16_ele,dipole_ent_s(m)));
    R26=abs(interp1(s_ele,R26_ele,dipole_exit_s(m+1))-interp1(s_ele,R26_ele,dipole_ent_s(m)));
    R56=abs(interp1(s_ele,R56_ele,dipole_exit_s(m+1))-interp1(s_ele,R56_ele,dipole_ent_s(m)));
    %R56_m=abs(interp1(s_ele,R56_ele,s_m(m))-interp1(s_ele,R56_ele,s_m(m+1)));
    Psi_fi=abs(psix_f(m+1)-psix_i(m));
    %Psi_fi=abs(psix_m(m)-psix_m(m+1));
    
    fprintf('from %d-th to %d-th dipole,R16=%f m,R26=%f,R56_m=%f m,Psix=%f rad \n',m,m+1,R16,R26,R56,Psi_fi);
    
    % scan betax to obtain mesh plot R56s1tos2 vs. Psi vs. betax
    for n=1:1:scan_num        
        R56s1tos2_01(n,:)=R56s1tos2_Plot3D(alphax_i(m),alphax_f(m+1),betax_i(m),betax_scan(n),Psi_scan,s1,s2,Lb,rho1,rho2,R16,R26,R56);
    end
    
    if (interp2(Psi_scan/pi,betax_scan,abs(R56s1tos2_01),Psi_fi/pi,betax_i(m))>=0.8*max(max(abs(R56s1tos2_01))))
        flag01=flag01+1;
    end
    
    % scan alphax to obtain mesh plot R56s1tos2 vs. Psi vs. alphax
    for n=1:1:scan_num        
        R56s1tos2_02(n,:)=R56s1tos2_Plot3D(alphax_i(m),alphax_scan(n),betax_i(m),betax_f(m+1),Psi_scan,s1,s2,Lb,rho1,rho2,R16,R26,R56);
    end
    
    if (interp2(Psi_scan/pi,alphax_scan,abs(R56s1tos2_02),Psi_fi/pi,alphax_i(m))>=0.8*max(max(abs(R56s1tos2_02))))
        flag02=flag02+1;
    end
    
    bottom=min(min(min(abs(R56s1tos2_01))),min(min(abs(R56s1tos2_02))));
    top=max(max(max(R56s1tos2_01)),max(max(R56s1tos2_02)));
    
    if (iplot==1)
    figure(m);
    h1=subplot(1,2,1); set(gca,'FontSize',40,'linewidth',5); plot3(Psi_fi/pi,betax_i(m),max(max(abs(R56s1tos2_01))),'kx','markersize',30,'linewidth',5); grid off; hold on;
    %subplot(1,2,1); set(gca,'FontSize',40,'linewidth',5); scatter3(Psi_fi,betax_i(m),max(max(abs(R56s1tos2_01))),'filled','MarkerFaceColor','black');
    h2=subplot(1,2,2); set(gca,'FontSize',40,'linewidth',5); plot3(Psi_fi/pi,alphax_i(m),max(max(abs(R56s1tos2_02))),'kx','markersize',30,'linewidth',5); grid off; hold on;
    
    figure(m);
    h1=subplot(1,2,1); set(gca,'FontSize',40,'linewidth',5); surf(Psi_scan/pi,betax_scan,abs(R56s1tos2_01)); xlabel('\Psi_{21} (\pi)'); ylabel('\beta_1 (m)'); xlim([0 2]); shading interp; hold on; view(90,-90); alpha(0.8); h=colorbar; set(h,'fontsize',30); axis('tight'); %caxis manual; caxis([bottom top]);
    h2=subplot(1,2,2); set(gca,'FontSize',40,'linewidth',5); surf(Psi_scan/pi,alphax_scan,abs(R56s1tos2_02)); xlabel('\Psi_{21} (\pi)'); ylabel('\alpha_1 (m)'); xlim([0 2]); shading interp; hold on; view(90,-90); alpha(0.8); h=colorbar; set(h,'fontsize',30); axis('tight'); %caxis manual; caxis([bottom top]); colorbar;
    end
end

fprintf('Dangerous ratio %f%%, %f%%...\n',flag01/r*100,flag02/r*100);
fprintf('flag01=%d, flag02=%d, num_r=%d \n',flag01,flag02,r);

%}
%=========================================================================%
%
[r,~]=size(s_m);
scan_num=40;

flag01=0; flag02=0; num_r=0;
for m=1:1:(r-1)
    for p=m+1:1:r
    num_r=num_r+1;
    alphax_scan=linspace(0.1*alphax_i(m),10*alphax_i(m),scan_num);
    betax_scan=linspace(0.1*betax_i(m),5*betax_i(m),scan_num);
    Psi_scan=linspace(0,2*pi,20);
    s1=0.5*abs(dipole_ent_s(m)-dipole_exit_s(m));
    s2=0.5*abs(dipole_ent_s(p)-dipole_exit_s(p));
    Lb=2*s1;
    rho1=dipole_mat(2*m,2);
    rho2=dipole_mat(2*p,2);
    
    if (rho1*rho2<0)
        s2=-s2;
        fprintf('%d-th to %d-th dipole reverse bend detected\n',m,p);
    end
    
    rho1=abs(rho1);
    rho2=abs(rho2);
    
    R16=abs(interp1(s_ele,R16_ele,dipole_exit_s(p))-interp1(s_ele,R16_ele,dipole_ent_s(m)));
    R26=abs(interp1(s_ele,R26_ele,dipole_exit_s(p))-interp1(s_ele,R26_ele,dipole_ent_s(m)));
    R56=abs(interp1(s_ele,R56_ele,dipole_exit_s(p))-interp1(s_ele,R56_ele,dipole_ent_s(m)));
    Psi_fi=abs(psix_f(p)-psix_i(m));
    
    fprintf('from %d-th to %d-th dipole,R16=%f m,R26=%f,R56=%f m,Psix=%f rad \n',m,p,R16,R26,R56,Psi_fi);
    
    % scan betax to obtain mesh plot R56s1tos2 vs. Psi vs. betax
    for n=1:1:scan_num        
        R56s1tos2_01(n,:)=R56s1tos2_Plot3D(alphax_i(m),alphax_f(p),betax_i(m),betax_scan(n),Psi_scan,s1,s2,Lb,rho1,rho2,R16,R26,R56);
    end
    
    if (interp2(Psi_scan/pi,betax_scan,abs(R56s1tos2_01),Psi_fi/pi,betax_i(m))>=0.8*max(max(abs(R56s1tos2_01))))
        flag01=flag01+1;
    end
    
    % scan alphax to obtain mesh plot R56s1tos2 vs. Psi vs. alphax
    for n=1:1:scan_num        
        R56s1tos2_02(n,:)=R56s1tos2_Plot3D(alphax_i(m),alphax_scan(n),betax_i(m),betax_f(p),Psi_scan,s1,s2,Lb,rho1,rho2,R16,R26,R56);
    end
    
    if (interp2(Psi_scan/pi,alphax_scan,abs(R56s1tos2_02),Psi_fi/pi,alphax_i(m))>=0.8*max(max(abs(R56s1tos2_02))))
        flag02=flag02+1;
    end
    
    bottom=min(min(min(abs(R56s1tos2_01))),min(min(abs(R56s1tos2_02))));
    top=max(max(max(R56s1tos2_01)),max(max(R56s1tos2_02)));
    
    if (iplot==1)
    figure(m*100+p);
    h1=subplot(1,2,1); set(gca,'FontSize',40,'linewidth',5); plot3(mod(Psi_fi,2*pi)/pi,betax_i(m),10*max(max(abs(R56s1tos2_01))),'kx','markersize',30,'linewidth',5); grid off; hold on;
    %subplot(1,2,1); set(gca,'FontSize',40,'linewidth',5); scatter3(Psi_fi,betax_i(m),max(max(abs(R56s1tos2_01))),'filled','MarkerFaceColor','black');
    h2=subplot(1,2,2); set(gca,'FontSize',40,'linewidth',5); plot3(mod(Psi_fi,2*pi)/pi,alphax_i(m),10*max(max(abs(R56s1tos2_02))),'kx','markersize',30,'linewidth',5); grid off; hold on;
    
    figure(m*100+p);
    h1=subplot(1,2,1); set(gca,'FontSize',40,'linewidth',5); surf(Psi_scan/pi,betax_scan,abs(R56s1tos2_01)); xlabel('\Psi_{21} (\pi)'); ylabel('\beta_1 (m)'); xlim([0 2]); shading interp; hold on; view(90,-90); alpha(0.8); h=colorbar; set(h,'fontsize',30); axis('tight'); %caxis manual; caxis([bottom top]);
    h2=subplot(1,2,2); set(gca,'FontSize',40,'linewidth',5); surf(Psi_scan/pi,alphax_scan,abs(R56s1tos2_02)); xlabel('\Psi_{21} (\pi)'); ylabel('\alpha_1 (m)'); xlim([0 2]); shading interp; hold on; view(90,-90); alpha(0.8); h=colorbar; set(h,'fontsize',30); axis('tight'); %caxis manual; caxis([bottom top]); colorbar;
    end
    
    end
end

fprintf('Dangerous ratio %f%%, %f%%...\n',flag01/num_r*100,flag02/num_r*100);
fprintf('flag01=%d, flag02=%d, num_r=%d \n',flag01,flag02,num_r);

%}


