if (scan_num==1)
fname='gain_func.sdds';

s_SI_unit=s/100;
d2d=abs(gkd_mat);
d2e=abs(ekd_mat);

fid=fopen(fname,'w');
fprintf(fid,'SDDS1\n');
fprintf(fid,'&parameter name=Description, type=string &end\n');
fprintf(fid,'&parameter name=I_b, type=double, description="Initial peak current in Amp" &end\n');
fprintf(fid,'&parameter name=C_ele, type=double, description="Total compression factor" &end\n');
fprintf(fid,'&parameter name=mesh_num, type=double, description="s-mesh number" &end\n');
fprintf(fid,'&parameter name=lambda, type=double, description="Initial mod. wavelength in um" &end\n');
fprintf(fid,'&parameter name=chirp, type=double, description="Initial chirp parameter in m^-1 " &end\n');
fprintf(fid,'&parameter name=iCSR_ss, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iCSR_tr, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iCSR_drift, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iLSC, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=ilinac, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=issCSRpp, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iIBS, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=itransLD, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iDz, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iDzz, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=DO_X, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=DO_Y, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=DO_Z, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=Derbenev, type=double, description=" " &end\n');
fprintf(fid,'&column name=s, units=m, type=double, description="Lattice s" &end\n');
fprintf(fid,'&column name=G11, units=none, type=double, description="Density-to-density gain" &end\n');
fprintf(fid,'&column name=G31, units=none, type=double, description="Density-to-energy gain" &end\n');
fprintf(fid,'&data mode=ascii, &end\n');
fprintf(fid,'X\n');
fprintf(fid,'%12.8e\n',I_b);
fprintf(fid,'%12.8e\n',C_ele(end));
fprintf(fid,'%12.8e\n',mesh_num);
fprintf(fid,'%12.8e\n',lambda_start01*1e4);
fprintf(fid,'%12.8e\n',chirp);
fprintf(fid,'%d\n',iCSR_ss);
fprintf(fid,'%d\n',iCSR_tr);
fprintf(fid,'%d\n',iCSR_drift);
fprintf(fid,'%d\n',iLSC);
fprintf(fid,'%d\n',ilinac);
fprintf(fid,'%d\n',issCSRpp);
fprintf(fid,'%d\n',iIBS);
fprintf(fid,'%d\n',itransLD);
fprintf(fid,'%d\n',iDz);
fprintf(fid,'%d\n',iDzz);
fprintf(fid,'%d\n',DO_X);
fprintf(fid,'%d\n',DO_Y);
fprintf(fid,'%d\n',DO_Z);
fprintf(fid,'%d\n',Derbenev);

fprintf(fid,'! page number 1\n');
fprintf(fid,'               %5.0f\n',length(s_SI_unit)-1);

for j=1:1:length(s_SI_unit)-1
	fprintf(fid,'%12.8e  %12.8e  %12.8e \n',s_SI_unit(j),d2d(j),d2e(j));
end

fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (scan_num>=2)
fname='gain_spec.sdds';

lambda_array_um=lambda_array*10^4;

fid=fopen(fname,'w');
fprintf(fid,'SDDS1\n');
fprintf(fid,'&parameter name=Description, type=string &end\n');
fprintf(fid,'&parameter name=I_b, type=double, description="Initial peak current in Amp" &end\n');
fprintf(fid,'&parameter name=C_ele, type=double, description="Total compression factor" &end\n');
fprintf(fid,'&parameter name=mesh_num, type=double, description="s-mesh number" &end\n');
fprintf(fid,'&parameter name=chirp, type=double, description="Initial chirp parameter in m^-1 " &end\n');
fprintf(fid,'&parameter name=iCSR_ss, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iCSR_tr, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iCSR_drift, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iLSC, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=ilinac, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=issCSRpp, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iIBS, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=itransLD, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iDz, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=iDzz, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=DO_X, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=DO_Y, type=double, description=" " &end\n');
fprintf(fid,'&parameter name=DO_Z, type=double, description=" " &end\n');
fprintf(fid,'&column name=lambda, units=um, type=double, description="Initial mod. wavelength in um" &end\n');
fprintf(fid,'&column name=Gf_11, units=none, type=double, description="Density-to-density gain" &end\n');
fprintf(fid,'&column name=Gf_31, units=none, type=double, description="Density-to-energy gain" &end\n');
fprintf(fid,'&data mode=ascii, &end\n');
fprintf(fid,'X\n');
fprintf(fid,'%12.8e\n',I_b);
fprintf(fid,'%12.8e\n',C_ele(end));
fprintf(fid,'%12.8e\n',mesh_num);
fprintf(fid,'%12.8e\n',chirp);
fprintf(fid,'%d\n',iCSR_ss);
fprintf(fid,'%d\n',iCSR_tr);
fprintf(fid,'%d\n',iCSR_drift);
fprintf(fid,'%d\n',iLSC);
fprintf(fid,'%d\n',ilinac);
fprintf(fid,'%d\n',issCSRpp);
fprintf(fid,'%d\n',iIBS);
fprintf(fid,'%d\n',itransLD);
fprintf(fid,'%d\n',iDz);
fprintf(fid,'%d\n',iDzz);
fprintf(fid,'%d\n',DO_X);
fprintf(fid,'%d\n',DO_Y);
fprintf(fid,'%d\n',DO_Z);

fprintf(fid,'! page number 1\n');
fprintf(fid,'               %5.0f\n',length(lambda_array));

for j=1:1:length(lambda_array)
	fprintf(fid,'%12.8e  %12.8e  %12.8e \n',lambda_array_um(j),Gf_11(j),Gf_31(j));
end

fclose(fid);
end
