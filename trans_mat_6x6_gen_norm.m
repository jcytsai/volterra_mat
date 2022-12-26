% ========================================================================%
% This program intends to derive 6x6 transport matrix from s to s' by
% tracking 36 individual/independent particles for ramping-energy beamline.
% Theoretical formulation can be found in the appendix of M. Venturini, 
% PRST-AB 10, 104401 (2007), based on re-defined canonical variables.
% Before running this program, one should use ELEGANT to track particles
% along the beamline while add individual WATCH points.
%
% Input: num_watch_pts, w001.o w002.o ...
%                       w001.s w002.s ...
%        pCentral as a function of s
%
% Output: 6x6 transport matrix as a function of s
%
% Program written by Cheng-Ying Tsai, jcytsai@vt.edu
% ========================================================================%
%clear; clc;
format long

c_speed=2.99792458e8;
num_watch_pts=TNWP;
num_particles=TNMP;
ratio=num_particles/36;
if (rem(ratio,1)~=0)
    fprintf('error: number of particles must be integer multiple of 36...\n');
    pause;
end
num_groups=num_particles/6;
if (rem(num_groups,1)~=0)
    fprintf('error: number of sub-matrices or groups must be an integer...\n');
    pause;
end
neg_sign=1;

% read 6D 4-particle distributions (w1,w2,...)
% format (x,xp,y,yp,t,p=beta*egamma)
% row: particle index (1-36)
% col: individual coord x,xp,y,yp,t,p
s_vec=zeros(1,num_watch_pts);
Z_map=zeros(num_particles,6,num_watch_pts);
for m=1:1:num_watch_pts
    
    filename01=sprintf('w%04ld.o', m);
    filename02=sprintf('w%04ld.s', m);
    
    % read 6-D coordinates (Z) at various locations as different pages (row x col x page)
    tmp01=load(filename01);
    
    % assign 2D array tmp01 to 3D array Z_map
    [r,c]=size(tmp01);
    for p=1:1:r
        for q=1:1:c
            Z_map(p,q,m)=tmp01(p,q);
        end
    end
    
    % read s-location for each WATCH point, in m
    s_vec(m)=load(filename02);
    
end


% read pCentral information
tmp00=load('pcentral_function.o');
pCentral=tmp00(:,2);
ECentral_tmp=0.511*sqrt(1+pCentral.^2);
ECentral=interp1(tmp00(:,1),ECentral_tmp,s_vec);
E0=ECentral(1);



% convert (x,xp,y,yp,t,p=beta*egamma) to (x,xp,y,yp,z,delta)
% z=(-)ct >0 for bunch head
% delta=(E-E0)/E0~(betagamma-mean_betagamma)/mean_betagamma
[r,c,pg]=size(Z_map);
mean_z=zeros(1,pg);
mean_betagamma=zeros(1,pg);
for m=1:1:pg
    %mean_x(m)=mean(Z_map(:,1,m));
    %mean_xp(m)=mean(Z_map(:,2,m));
    %mean_y(m)=mean(Z_map(:,3,m));
    %mean_yp(m)=mean(Z_map(:,4,m));
    mean_z(m)=neg_sign*mean(Z_map(:,5,m))*c_speed;
    mean_betagamma(m)=mean(Z_map(:,6,m));
    
    Z_map(:,1,m)=Z_map(:,1,m)*sqrt(ECentral(m)/E0);
    Z_map(:,2,m)=Z_map(:,2,m)*sqrt(ECentral(m)/E0);
    Z_map(:,3,m)=Z_map(:,3,m)*sqrt(ECentral(m)/E0);
    Z_map(:,4,m)=Z_map(:,4,m)*sqrt(ECentral(m)/E0);
    %Z_map(:,1,m)=(Z_map(:,1,m)-mean_x(m))*sqrt(ECentral(m)/E0);
    %Z_map(:,2,m)=(Z_map(:,2,m)-mean_xp(m))*sqrt(ECentral(m)/E0);
    %Z_map(:,3,m)=(Z_map(:,3,m)-mean_y(m))*sqrt(ECentral(m)/E0);
    %Z_map(:,4,m)=(Z_map(:,4,m)-mean_yp(m))*sqrt(ECentral(m)/E0);
    Z_map(:,5,m)=neg_sign*Z_map(:,5,m)*c_speed-mean_z(m);
    Z_map(:,6,m)=(Z_map(:,6,m)-mean_betagamma(m))/mean_betagamma(m);
    Z_map(:,6,m)=(Z_map(:,6,m)+1)-ECentral(m)/E0;
    
end


% prepare Z_i matrix
Zi_sub=squeeze(Z_map(:,:,1));
%Zi_mat=kron(eye(6),Zi_sub);
Zi_mat=zeros(num_particles,36);
for gp=1:1:num_groups
    if (rem(gp*6,36)==0) % when gp=6,12,18...
        tmp03=36;
        Zi_mat(gp*6-5:gp*6,tmp03-5:tmp03)=Zi_sub(gp*6-5:gp*6,1:6);
    else
        Zi_mat(gp*6-5:gp*6,rem(gp*6,36)-5:rem(gp*6,36))=Zi_sub(gp*6-5:gp*6,1:6);
    end
end



% prepare Z_f vector for each s
for m=1:1:pg
    tmp02=squeeze(Z_map(:,:,m));
    
    for gp=1:1:num_groups
        if (rem(gp,6)==0)
            tmp04=6;
            Zf_s(gp*6-5:gp*6,m)=tmp02(gp*6-5:gp*6,tmp04);
        else
            Zf_s(gp*6-5:gp*6,m)=tmp02(gp*6-5:gp*6,rem(gp,6));
        end
    end

end


% calculate R_vec in vector form & re-shape R_vec to R_mat
for m=1:1:num_watch_pts
    R_vec(:,m)=Zi_mat\Zf_s(:,m);
    %R_vec(:,m)=inv(Zi_mat)*Zf_s(:,m);
    tmp=reshape(R_vec(:,m),6,6);
    A(:,:,m)=tmp';
end

% plot results
figure(1);
subplot(2,2,1); plot(s_vec,squeeze(A(1,6,:)),'ro-'); xlabel('s (m)'); ylabel('R_{16} (m)'); title('red-tracking; blue-matrix output'); hold on;
subplot(2,2,2); plot(s_vec,squeeze(A(5,1,:)),'ro-'); xlabel('s (m)'); ylabel('R_{51}'); hold on;
subplot(2,2,3); plot(s_vec,squeeze(A(5,2,:)),'ro-'); xlabel('s (m)'); ylabel('R_{52} (m)'); hold on;
subplot(2,2,4); plot(s_vec,squeeze(A(5,6,:)),'ro-'); xlabel('s (m)'); ylabel('R_{56} (m)'); hold on;

transport_ele=load('lattice_transport_functions_ELEGANT_original_matrix_output.o');
figure(1);
subplot(2,2,1); plot(transport_ele(:,1),transport_ele(:,7),'bo-'); hold on;
subplot(2,2,2); plot(transport_ele(:,1),transport_ele(:,26),'bo-'); hold on;
subplot(2,2,3); plot(transport_ele(:,1),transport_ele(:,27),'bo-'); hold on;
subplot(2,2,4); plot(transport_ele(:,1),transport_ele(:,31),'bo-'); hold on;

%pause;

% write results to file
% format (s,R11,R12,...R65,R66) with matrix size (num_watch_pts)x(1+36)
transport(:,1)=s_vec';
for m=1:1:num_watch_pts
    transport(m,2:37)=reshape(A(:,:,m)',1,36);
end

fileID1=fopen('lattice_transport_functions_ELEGANT.o','w');
fprintf(fileID1,'%10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f\n',transport');
fclose(fileID1);
