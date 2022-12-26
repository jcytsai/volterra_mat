function func=R56s1tos2_Plot3D(alphax_i,alphax_f,betax_i,betax_f,Psi_vec,s1,s2,Lb,rho1,rho2,R16,R26,R56)

% Note: 
% quantities with subscript i means entrance to 1st dipole
% qunatities with subscript f means exit of 2nd dipole

%tmp01=(Lb-s2)/rho2;
tmp01=(s2-Lb)/rho2;
tmp02=s1/rho1;
tmp03=sqrt(betax_i*betax_f);
tmp04=sqrt(betax_i/betax_f);
tmp05=sqrt(betax_f/betax_i);
tmp06=cos(Psi_vec)-alphax_f*sin(Psi_vec);
tmp07=cos(Psi_vec)+alphax_i*sin(Psi_vec);
tmp08=(alphax_f-alphax_i)*cos(Psi_vec)+(1+alphax_i*alphax_f)*sin(Psi_vec);

%tmp09=R56-Lb+s1+s2;
tmp09=R56+Lb+s1-s2;
tmp10=R26*rho2*(cos(tmp01)-1)-rho1*sin(tmp02)-R16*sin(tmp01)+rho2*sin(tmp01);
tmp11=-R26*tmp03*sin(Psi_vec)-tmp03*sin(tmp01)*sin(Psi_vec);
tmp12=R16*tmp04*tmp06+tmp04*rho2*(cos(tmp01)-1)*tmp06;
tmp13=sin(tmp02)*(tmp11+tmp12);
tmp14=-R26*tmp05*tmp07-tmp05*sin(tmp01)*tmp07;
tmp15=-R16/tmp03*tmp08-rho2*(cos(tmp01)-1)/tmp03*tmp08;
tmp16=rho1*(1-cos(tmp02))*(tmp14+tmp15);

func=tmp09+tmp10-tmp13+tmp16; % row vector of the same size with Psi_vec
