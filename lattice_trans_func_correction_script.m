%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This subroutine loads the lattice transport functions output from       %
% ELEGANT, read, check if duplicate s coordinates are present, delete     %
% them, and write an output corrected file with file name                 %
% lattice_transport_functions.o                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;
clearvars -EXCEPT shortest;
format long;

%============ load the ELEGANT lattice transport functions ===============%
filename='lattice_transport_functions_ELEGANT.o';
% format (s,R11,R12,...R16,R21,...,R26,...,R51,R52,R53,R54,R55,R56,...,R65,R66)
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
    fprintf('writing the corrected lattice transport functions... \n');
else
    fprintf('there is no duplicate s coordinate in the file... \n');
end

fprintf('generating the file lattice_transport_functions.o...');
fileID1 = fopen('lattice_transport_functions.o','w');
fprintf(fileID1,'%10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f %10.15f \n',transport_new');
fclose(fileID1);

fprintf('done...\n');
%=========================================================================%