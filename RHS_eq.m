function func=RHS_eq(s,Y,E0)

global alpha_RF lambda_RF k_RF omega_RF psi0;

c_speed=2.99792458e8;                 % unit: m/sec

tmp_t=Y(1);
tmp_egamma=Y(2);

tmp_M11=Y(3);
tmp_M12=Y(4);
tmp_M13=Y(5);
tmp_M14=Y(6);
tmp_M15=Y(7);
tmp_M16=Y(8);

tmp_M21=Y(9);
tmp_M22=Y(10);
tmp_M23=Y(11);
tmp_M24=Y(12);
tmp_M25=Y(13);
tmp_M26=Y(14);

tmp_M31=Y(15);
tmp_M32=Y(16);
tmp_M33=Y(17);
tmp_M34=Y(18);
tmp_M35=Y(19);
tmp_M36=Y(20);

tmp_M41=Y(21);
tmp_M42=Y(22);
tmp_M43=Y(23);
tmp_M44=Y(24);
tmp_M45=Y(25);
tmp_M46=Y(26);

tmp_M51=Y(27);
tmp_M52=Y(28);
tmp_M53=Y(29);
tmp_M54=Y(30);
tmp_M55=Y(31);
tmp_M56=Y(32);

tmp_M61=Y(33);
tmp_M62=Y(34);
tmp_M63=Y(35);
tmp_M64=Y(36);
tmp_M65=Y(37);
tmp_M66=Y(38);

tmp_s=Y(39);

beta_rel=sqrt(1-1/tmp_egamma^2);
%k_RF=2*pi/beta_rel/lambda_RF;
k_RF=2*pi/lambda_RF;
alpha_RF=E0/0.511/k_RF;

tmp_psi=k_RF*tmp_s-omega_RF*tmp_t+psi0;
%tmp_psi=k_RF*s-omega_RF*tmp_t+psi0;
%tmp_psi=psi0;

A12=1;
%A21=-alpha_RF*k_RF^2/(2*tmp_egamma)*cos(tmp_psi);
A21=-alpha_RF*k_RF^2/2*cos(tmp_psi);
A34=1;
A43=A21;
%A56=(tmp_egamma^2-1)^(-3/2);
A56=1/(beta_rel^2*tmp_egamma^3);
A65=-alpha_RF*k_RF^2*cos(tmp_psi);


tmp01=tmp_egamma/c_speed/sqrt(tmp_egamma^2-1);
tmp02=-alpha_RF*k_RF*sin(tmp_psi);
tmp03=tmp_M21;
tmp04=tmp_M22;
tmp05=tmp_M23;
tmp06=tmp_M24;
tmp07=tmp_M25;
tmp08=tmp_M26;
tmp09=A21*tmp_M11;
tmp10=A21*tmp_M12;
tmp11=A21*tmp_M13;
tmp12=A21*tmp_M14;
tmp13=A21*tmp_M15;
tmp14=A21*tmp_M16;
tmp15=tmp_M41;
tmp16=tmp_M42;
tmp17=tmp_M43;
tmp18=tmp_M44;
tmp19=tmp_M45;
tmp20=tmp_M46;
tmp21=A43*tmp_M31;
tmp22=A43*tmp_M32;
tmp23=A43*tmp_M33;
tmp24=A43*tmp_M34;
tmp25=A43*tmp_M35;
tmp26=A43*tmp_M36;
tmp27=A56*tmp_M61;
tmp28=A56*tmp_M62;
tmp29=A56*tmp_M63;
tmp30=A56*tmp_M64;
tmp31=A56*tmp_M65;
tmp32=A56*tmp_M66;
tmp33=A65*tmp_M51;
tmp34=A65*tmp_M52;
tmp35=A65*tmp_M53;
tmp36=A65*tmp_M54;
tmp37=A65*tmp_M55;
tmp38=A65*tmp_M56;
tmp39=1;


func=[tmp01;tmp02;...
      tmp03;tmp04;tmp05;tmp06;tmp07;tmp08;...
      tmp09;tmp10;tmp11;tmp12;tmp13;tmp14;...
      tmp15;tmp16;tmp17;tmp18;tmp19;tmp20;...
      tmp21;tmp22;tmp23;tmp24;tmp25;tmp26;...
      tmp27;tmp28;tmp29;tmp30;tmp31;tmp32;...
      tmp33;tmp34;tmp35;tmp36;tmp37;tmp38;...
      tmp39];





