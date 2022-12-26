function [func1,func2]=read_tot_col_vec_2x2(nrow,ncol)

file_gk=fopen('gk_mat_01.bin','r');
file_ek=fopen('ek_mat_01.bin','r');

tmp01=fread(file_gk,[nrow, ncol],'double'); tmp01=tmp01';
tmp02=fread(file_ek,[nrow, ncol],'double'); tmp02=tmp02';

fclose(file_gk);
fclose(file_ek);

gk_vec=complex(tmp01(:,1),tmp01(:,2));
ek_vec=complex(tmp02(:,1),tmp02(:,2));

func1=gk_vec;
func2=ek_vec;
