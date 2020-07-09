%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subroutine loads the pCentral function output from ELEGANT (*.cen),%
% read, check if duplicate s coordinates are present, delete them, and    %
% write an output corrected file with file name pcentral_function.o       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;
clearvars -EXCEPT shortest;
format long;

%============ load the ELEGANT lattice pcentral functions ===============%
filename='pcentral_function_ELEGANT.o';
% format (s,pCentral=beta*gamma)
delimiterIn=' '; headerlinesIn=0;
pcentral=importdata(filename,delimiterIn,headerlinesIn);

[r,~]=size(pcentral);
%=========================================================================%



%============ check if duplicate s coordinates are present ===============%
pcentral_new(1,:)=pcentral(1,:); n=2; q=0;
for m=2:1:r
    if (pcentral(m-1,1)==pcentral(m,1))
        %fprintf('duplicate s = %f m found in %d-th line... \n',pcentral(m,1),m);
        q=q+1;
    else
        pcentral_new(n,:)=pcentral(m,:); n=n+1;
    end
end
%=========================================================================%


%============ summary and output corrected lattice file ==================%
fprintf('total number of %d duplicate s coordinates are detected... \n',q);

if (q ~= 0)
    fprintf('writing the corrected pcentral data... \n');
else
    fprintf('there is no duplicate s coordinate in the file... \n');
end

fprintf('generating the file pcentral_function.o...');
fileID1 = fopen('pcentral_function.o','w');
fprintf(fileID1,'%10.15f %10.15f \n',pcentral_new');
fclose(fileID1);
fprintf('done...\n');
%=========================================================================%