%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subroutine loads the lattice element geometrical coordinates from  %
% ELEGANT, read, select dipole elements, and write an output file with    %
% file name lattice_dipole_info.o                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
format long;

%============ load the ELEGANT element geometrical information ===========%
filename='elements_geom.o';
% format (ElementType,s,theta,phi)
delimiterIn='\t'; headerlinesIn=0;
dipole_info=importdata(filename,delimiterIn,headerlinesIn);

[r,~]=size(dipole_info.data);
%=========================================================================%


%============ select dipole elements and read dipole information =========%
n=0; p=0; q=0; w=0;
for m=1:1:r
    
    if ((strcmp(dipole_info.textdata(m),'CSRCSBEND')==1)||...
        (strcmp(dipole_info.textdata(m),'CSBEND')==1)||...
        (strcmp(dipole_info.textdata(m),'KSBEND')==1)||...
        (strcmp(dipole_info.textdata(m),'NISBEND')==1)||...
        (strcmp(dipole_info.textdata(m),'RBEN')==1)||...
        (strcmp(dipole_info.textdata(m),'SBEN')==1)||...
        (strcmp(dipole_info.textdata(m),'TUBEND')==1))
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

dipole_ent_s=dipole_info.data(bend_ent_row_index,1);
dipole_exit_s=dipole_info.data(bend_ext_row_index,1);
dipole_beg_theta=dipole_info.data(bend_ent_row_index,2); % horizontal bend
dipole_end_theta=dipole_info.data(bend_ext_row_index,2);
dipole_beg_phi=dipole_info.data(bend_ent_row_index,3);   % vertical bend
dipole_end_phi=dipole_info.data(bend_ext_row_index,3);

%============ format conversion ==========================================%
dipole_s(1:2:2*n-1)=dipole_ent_s;
dipole_s(2:2:2*n)=dipole_exit_s;

dipole_theta(1:2:2*n-1)=dipole_beg_theta;
dipole_theta(2:2:2*n)=dipole_end_theta;

dipole_phi(1:2:2*n-1)=dipole_beg_phi;
dipole_phi(2:2:2*n)=dipole_end_phi;

for m=1:1:n
    dipole_rhox(2*m-1)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_theta(2*m))-(dipole_theta(2*m-1))); 
    if (((dipole_theta(2*m))-(dipole_theta(2*m-1))) > pi) dipole_rhox(2*m-1)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_theta(2*m))-(dipole_theta(2*m-1))-2*pi); end
    if (dipole_rhox(2*m-1)>1e6) dipole_rhox(2*m-1)=1e50; end
    
    dipole_rhox(2*m)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_theta(2*m))-(dipole_theta(2*m-1))); 
    if (((dipole_theta(2*m))-(dipole_theta(2*m-1))) > pi) dipole_rhox(2*m)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_theta(2*m))-(dipole_theta(2*m-1))-2*pi); end
    if (dipole_rhox(2*m)>1e6) dipole_rhox(2*m)=1e50; end
    
    dipole_rhoy(2*m-1)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_phi(2*m))-(dipole_phi(2*m-1))); 
    if (((dipole_phi(2*m))-(dipole_phi(2*m-1))) > pi) dipole_rhoy(2*m-1)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_phi(2*m))-(dipole_phi(2*m-1))-2*pi); end
    if (dipole_rhoy(2*m-1)>1e6) dipole_rhoy(2*m-1)=1e50; end
    
    dipole_rhoy(2*m)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_phi(2*m))-(dipole_phi(2*m-1))); 
    if (((dipole_phi(2*m))-(dipole_phi(2*m-1))) > pi) dipole_rhoy(2*m)=(dipole_s(2*m)-dipole_s(2*m-1))/((dipole_phi(2*m))-(dipole_phi(2*m-1))-2*pi); end
    if (dipole_rhoy(2*m)>1e6) dipole_rhoy(2*m)=1e50; end
    
end

dipole_mat=[dipole_s;dipole_rhox;dipole_rhoy]';
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

fileID1 = fopen('lattice_dipole_info.o','w');
fprintf(fileID1,'%5.5f %5.5f %5.5f \n',dipole_mat');
fclose(fileID1);

%=========================================================================%







