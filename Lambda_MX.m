function func=Lambda_MX(Lg_mod,wiggler_kL,wiggler_lambda,sig_x,beta_ave,emit_gx,sigma_delta)

a1=0.45; a2=0.57; a3=0.55; a4=1.6; a5=3;
a6=2; a7=0.35; a8=2.9; a9=2.4; a10=51;
a11=0.95; a12=3; a13=5.4; a14=0.7; a15=1.9;
a16=1140; a17=2.2; a18=2.9; a19=3.2;

% NOTE: gain length must be divided by sqrt(3) to be consistent with
% Kim, Huang, and Lindberg's book, p.179 (Eng. ed.)
eta_d=Lg_mod/(2*wiggler_kL*sig_x^2)/sqrt(3);
eta_e=2*Lg_mod*wiggler_kL*emit_gx/beta_ave/sqrt(3);
eta_r=4*pi*Lg_mod*sigma_delta/wiggler_lambda/sqrt(3);

tmp01=a1*eta_d^a2;
tmp02=a3*eta_e^a4;
tmp03=a5*eta_r^a6;
tmp04=a7*eta_e^a8*eta_r^a9;
tmp05=a10*eta_d^a11*eta_r^a12;
tmp06=a13*eta_d^a14*eta_e^a15;
tmp07=a16*eta_d^a17*eta_e^a18*eta_r^a19;

func=tmp01+tmp02+tmp03+tmp04+tmp05+tmp06+tmp07;