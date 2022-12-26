function [func1,func2,func3,func4,func5,func6,func7,func8,func9,func10,func11,func12,func13,func14,func15,func16,func17,func18,func19,func20,func21,func22,func23,func24,func25,func26,func27,func28,func29,func30,func31,func32,func33,func34,func35,func36]=read_tot_col_vec_36x36(nrow,ncol,file_index)

file_gkd    =fopen(sprintf('gkd_mat_%02d.bin',file_index),'r');
file_gkp    =fopen(sprintf('gkp_mat_%02d.bin',file_index),'r');
file_gkxz   =fopen(sprintf('gkxz_mat_%02d.bin',file_index),'r');
file_gkxpz  =fopen(sprintf('gkxpz_mat_%02d.bin',file_index),'r');
file_gkyz   =fopen(sprintf('gkyz_mat_%02d.bin',file_index),'r');
file_gkypz  =fopen(sprintf('gkypz_mat_%02d.bin',file_index),'r');

file_ekd    =fopen(sprintf('ekd_mat_%02d.bin',file_index),'r');
file_ekp    =fopen(sprintf('ekp_mat_%02d.bin',file_index),'r');
file_ekxz   =fopen(sprintf('ekxz_mat_%02d.bin',file_index),'r');
file_ekxpz  =fopen(sprintf('ekxpz_mat_%02d.bin',file_index),'r');
file_ekyz   =fopen(sprintf('ekyz_mat_%02d.bin',file_index),'r');
file_ekypz  =fopen(sprintf('ekypz_mat_%02d.bin',file_index),'r');

file_axkd   =fopen(sprintf('axkd_mat_%02d.bin',file_index),'r');
file_axkp   =fopen(sprintf('axkp_mat_%02d.bin',file_index),'r');
file_axkxz  =fopen(sprintf('axkxz_mat_%02d.bin',file_index),'r');
file_axkxpz =fopen(sprintf('axkxpz_mat_%02d.bin',file_index),'r');
file_axkyz  =fopen(sprintf('axkyz_mat_%02d.bin',file_index),'r');
file_axkypz =fopen(sprintf('axkypz_mat_%02d.bin',file_index),'r');

file_axpkd  =fopen(sprintf('axpkd_mat_%02d.bin',file_index),'r');
file_axpkp  =fopen(sprintf('axpkp_mat_%02d.bin',file_index),'r');
file_axpkxz =fopen(sprintf('axpkxz_mat_%02d.bin',file_index),'r');
file_axpkxpz=fopen(sprintf('axpkxpz_mat_%02d.bin',file_index),'r');
file_axpkyz =fopen(sprintf('axpkxz_mat_%02d.bin',file_index),'r');
file_axpkypz=fopen(sprintf('axpkxpz_mat_%02d.bin',file_index),'r');

file_aykd   =fopen(sprintf('aykd_mat_%02d.bin',file_index),'r');
file_aykp   =fopen(sprintf('aykp_mat_%02d.bin',file_index),'r');
file_aykxz  =fopen(sprintf('aykxz_mat_%02d.bin',file_index),'r');
file_aykxpz =fopen(sprintf('aykxpz_mat_%02d.bin',file_index),'r');
file_aykyz  =fopen(sprintf('aykyz_mat_%02d.bin',file_index),'r');
file_aykypz =fopen(sprintf('aykypz_mat_%02d.bin',file_index),'r');

file_aypkd  =fopen(sprintf('aypkd_mat_%02d.bin',file_index),'r');
file_aypkp  =fopen(sprintf('aypkp_mat_%02d.bin',file_index),'r');
file_aypkxz =fopen(sprintf('aypkxz_mat_%02d.bin',file_index),'r');
file_aypkxpz=fopen(sprintf('aypkxpz_mat_%02d.bin',file_index),'r');
file_aypkyz =fopen(sprintf('aypkyz_mat_%02d.bin',file_index),'r');
file_aypkypz=fopen(sprintf('aypkypz_mat_%02d.bin',file_index),'r');

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
tmp03=fread(file_gkxz,[nrow, ncol],'double'); tmp03=tmp03';
tmp04=fread(file_gkxpz,[nrow, ncol],'double'); tmp04=tmp04';
tmp05=fread(file_gkyz,[nrow, ncol],'double'); tmp05=tmp05';
tmp06=fread(file_gkypz,[nrow, ncol],'double'); tmp06=tmp06';

tmp07=fread(file_ekd,[nrow, ncol],'double'); tmp07=tmp07';
tmp08=fread(file_ekp,[nrow, ncol],'double'); tmp08=tmp08';
tmp09=fread(file_ekxz,[nrow, ncol],'double'); tmp09=tmp09';
tmp10=fread(file_ekxpz,[nrow, ncol],'double'); tmp10=tmp10';
tmp11=fread(file_ekyz,[nrow, ncol],'double'); tmp11=tmp11';
tmp12=fread(file_ekypz,[nrow, ncol],'double'); tmp12=tmp12';

tmp13=fread(file_axkd,[nrow, ncol],'double'); tmp13=tmp13';
tmp14=fread(file_axkp,[nrow, ncol],'double'); tmp14=tmp14';
tmp15=fread(file_axkxz,[nrow, ncol],'double'); tmp15=tmp15';
tmp16=fread(file_axkxpz,[nrow, ncol],'double'); tmp16=tmp16';
tmp17=fread(file_axkyz,[nrow, ncol],'double'); tmp17=tmp17';
tmp18=fread(file_axkypz,[nrow, ncol],'double'); tmp18=tmp18';

tmp19=fread(file_axpkd,[nrow, ncol],'double'); tmp19=tmp19';
tmp20=fread(file_axpkp,[nrow, ncol],'double'); tmp20=tmp20';
tmp21=fread(file_axpkxz,[nrow, ncol],'double'); tmp21=tmp21';
tmp22=fread(file_axpkxpz,[nrow, ncol],'double'); tmp22=tmp22';
tmp23=fread(file_axpkyz,[nrow, ncol],'double'); tmp23=tmp23';
tmp24=fread(file_axpkypz,[nrow, ncol],'double'); tmp24=tmp24';

tmp25=fread(file_aykd,[nrow, ncol],'double'); tmp25=tmp25';
tmp26=fread(file_aykp,[nrow, ncol],'double'); tmp26=tmp26';
tmp27=fread(file_aykxz,[nrow, ncol],'double'); tmp27=tmp27';
tmp28=fread(file_aykxpz,[nrow, ncol],'double'); tmp28=tmp28';
tmp29=fread(file_aykyz,[nrow, ncol],'double'); tmp29=tmp29';
tmp30=fread(file_aykypz,[nrow, ncol],'double'); tmp30=tmp30';

tmp31=fread(file_aypkd,[nrow, ncol],'double'); tmp31=tmp31';
tmp32=fread(file_aypkp,[nrow, ncol],'double'); tmp32=tmp32';
tmp33=fread(file_aypkxz,[nrow, ncol],'double'); tmp33=tmp33';
tmp34=fread(file_aypkxpz,[nrow, ncol],'double'); tmp34=tmp34';
tmp35=fread(file_aypkyz,[nrow, ncol],'double'); tmp35=tmp35';
tmp36=fread(file_aypkypz,[nrow, ncol],'double'); tmp36=tmp36';

fclose('all');

gkd_vec=complex(tmp01(:,1),tmp01(:,2));
gkp_vec=complex(tmp02(:,1),tmp02(:,2));
gkxz_vec=complex(tmp03(:,1),tmp03(:,2));
gkxpz_vec=complex(tmp04(:,1),tmp04(:,2));
gkyz_vec=complex(tmp05(:,1),tmp05(:,2));
gkypz_vec=complex(tmp06(:,1),tmp06(:,2));

ekd_vec=complex(tmp07(:,1),tmp07(:,2));
ekp_vec=complex(tmp08(:,1),tmp08(:,2));
ekxz_vec=complex(tmp09(:,1),tmp09(:,2));
ekxpz_vec=complex(tmp10(:,1),tmp10(:,2));
ekyz_vec=complex(tmp11(:,1),tmp11(:,2));
ekypz_vec=complex(tmp12(:,1),tmp12(:,2));

axkd_vec=complex(tmp13(:,1),tmp13(:,2));
axkp_vec=complex(tmp14(:,1),tmp14(:,2));
axkxz_vec=complex(tmp15(:,1),tmp15(:,2));
axkxpz_vec=complex(tmp16(:,1),tmp16(:,2));
axkyz_vec=complex(tmp17(:,1),tmp17(:,2));
axkypz_vec=complex(tmp18(:,1),tmp18(:,2));

axpkd_vec=complex(tmp19(:,1),tmp19(:,2));
axpkp_vec=complex(tmp20(:,1),tmp20(:,2));
axpkxz_vec=complex(tmp21(:,1),tmp21(:,2));
axpkxpz_vec=complex(tmp22(:,1),tmp22(:,2));
axpkyz_vec=complex(tmp23(:,1),tmp23(:,2));
axpkypz_vec=complex(tmp24(:,1),tmp24(:,2));

aykd_vec=complex(tmp25(:,1),tmp25(:,2));
aykp_vec=complex(tmp26(:,1),tmp26(:,2));
aykxz_vec=complex(tmp27(:,1),tmp27(:,2));
aykxpz_vec=complex(tmp28(:,1),tmp28(:,2));
aykyz_vec=complex(tmp29(:,1),tmp29(:,2));
aykypz_vec=complex(tmp30(:,1),tmp30(:,2));

aypkd_vec=complex(tmp31(:,1),tmp31(:,2));
aypkp_vec=complex(tmp32(:,1),tmp32(:,2));
aypkxz_vec=complex(tmp33(:,1),tmp33(:,2));
aypkxpz_vec=complex(tmp34(:,1),tmp34(:,2));
aypkyz_vec=complex(tmp35(:,1),tmp35(:,2));
aypkypz_vec=complex(tmp36(:,1),tmp36(:,2));

%tot_col_vec_old=[gkd_vec gkp_vec ekd_vec ekp_vec];
%tot_col_vec_old=[gkd_vec ekp_vec];
%func=tot_col_vec_old;
func1=gkd_vec;
func2=gkp_vec;
func3=gkxz_vec;
func4=gkxpz_vec;
func5=gkyz_vec;
func6=gkypz_vec;

func7=ekd_vec;
func8=ekp_vec;
func9=ekxz_vec;
func10=ekxpz_vec;
func11=ekyz_vec;
func12=ekypz_vec;

func13=axkd_vec;
func14=axkp_vec;
func15=axkxz_vec;
func16=axkxpz_vec;
func17=axkyz_vec;
func18=axkypz_vec;

func19=axpkd_vec;
func20=axpkp_vec;
func21=axpkxz_vec;
func22=axpkxpz_vec;
func23=axkyz_vec;
func24=axkypz_vec;

func25=aykd_vec;
func26=aykp_vec;
func27=aykxz_vec;
func28=aykxpz_vec;
func29=aykyz_vec;
func30=aykypz_vec;

func31=aypkd_vec;
func32=aypkp_vec;
func33=aypkxz_vec;
func34=aypkxpz_vec;
func35=aypkyz_vec;
func36=aypkypz_vec;

