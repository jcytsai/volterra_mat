%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subroutine loads the beam sigma matrix functions output from       %
% ELEGANT, read, check if duplicate s coordinates are present, delete     %
% them, and write an output corrected file with file name                 %
% beam_sigma_mat_functions.o                                              %
% 20200225: two more columns are added, Sxp & Syp                         % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
clearvars -EXCEPT shortest;
format long

%=============== load the ELEGANT sigma matrix functions =================%
filename='beam_sigma_mat_functions_ELEGANT.o';
% format (s,Sx,Sy,Ss,Sdelta,enx,eny,Sxp,Syp)
delimiterIn=' '; headerlinesIn=0;
transport=importdata(filename,delimiterIn,headerlinesIn);

[r,~]=size(transport);
%=========================================================================%



%============ check if duplicate s coordinates are present ===============%
transport_new(1,:)=transport(1,:); n=2; q=0;
for m=2:1:r
    if (transport(m-1,1)==transport(m,1))
        %fprintf('duplicate s = %f m found in %d-th line... \n',transport(m,1),m);
        q=q+1;
    else
        transport_new(n,:)=transport(m,:); n=n+1;
    end
end
%=========================================================================%


%============ summary and output corrected lattice file ==================%
fprintf('total number of %d duplicate s coordinates are detected... \n',q);

if (q ~= 0)
    fprintf('writing the corrected beam sigma matrix functions... \n');
else
    fprintf('there is no duplicate s coordinate in the file... \n');
end

fprintf('generating the file beam_sigma_mat_functions.o...');
fileID1 = fopen('beam_sigma_mat_functions.o','w');
fprintf(fileID1,'%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f\n',transport_new');
fclose(fileID1);
fprintf('done...\n');
%=========================================================================%


