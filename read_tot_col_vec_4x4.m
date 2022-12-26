function [func1,func2,func3,func4]=read_tot_col_vec_4x4(nrow,ncol,file_index)

file_gkd=fopen(sprintf('gkd_mat_%02d.bin',file_index),'r');
file_gkp=fopen(sprintf('gkp_mat_%02d.bin',file_index),'r');
file_ekd=fopen(sprintf('ekd_mat_%02d.bin',file_index),'r');
file_ekp=fopen(sprintf('ekp_mat_%02d.bin',file_index),'r');

% format:
% input [Re(1) Im(1) Re(2) Im(2) Re(3) Im(3)...]
% fread [Re(1) Re(2) Re(3)...]
%       [Im(1) Im(2) Im(3)...]
% transpose [Re(1) Im(1)]
%           [Re(2) Im(2)]
%           [Re(3) Im(3)]
%           [...   ...  ]

tmp01=fread(file_gkd,[nrow, ncol],'double'); tmp01=tmp01';
tmp02=fread(file_gkp,[nrow, ncol],'double'); tmp02=tmp02';
tmp03=fread(file_ekd,[nrow, ncol],'double'); tmp03=tmp03';
tmp04=fread(file_ekp,[nrow, ncol],'double'); tmp04=tmp04';

fclose(file_gkd);
fclose(file_gkp);
fclose(file_ekd);
fclose(file_ekp);

gkd_vec=complex(tmp01(:,1),tmp01(:,2));
gkp_vec=complex(tmp02(:,1),tmp02(:,2));
ekd_vec=complex(tmp03(:,1),tmp03(:,2));
ekp_vec=complex(tmp04(:,1),tmp04(:,2));

%tot_col_vec_old=[gkd_vec gkp_vec ekd_vec ekp_vec];
%tot_col_vec_old=[gkd_vec ekp_vec];
%func=tot_col_vec_old;
func1=gkd_vec;
func2=gkp_vec;
func3=ekd_vec;
func4=ekp_vec;
