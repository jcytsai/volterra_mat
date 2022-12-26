% this special program intends to extract stage gain coefficient d_m
% directly from G_f^(M), instead of using polynominal fitting technique.
% For more information, see Eqs.(10,11,14) of the long paper.

Alfven=17045;

coeff=const_A*I_b/(egamma*Alfven);
term=coeff.^(0:1:10);

gk_0=G0_k(end);             n1k=C_mat(end)*gk_0; G_00_c=n1k/(C_mat(end)*n_1k0);
gk_1=tmp01*transpose(G0_k); n1k=C_mat(end)*gk_1; G_01_c=n1k/(C_mat(end)*n_1k0);
gk_2=tmp02*transpose(G0_k); n1k=C_mat(end)*gk_2; G_02_c=n1k/(C_mat(end)*n_1k0);
gk_3=tmp03*transpose(G0_k); n1k=C_mat(end)*gk_3; G_03_c=n1k/(C_mat(end)*n_1k0);
gk_4=tmp04*transpose(G0_k); n1k=C_mat(end)*gk_4; G_04_c=n1k/(C_mat(end)*n_1k0);
gk_5=tmp05*transpose(G0_k); n1k=C_mat(end)*gk_5; G_05_c=n1k/(C_mat(end)*n_1k0);
gk_6=tmp06*transpose(G0_k); n1k=C_mat(end)*gk_6; G_06_c=n1k/(C_mat(end)*n_1k0);
gk_7=tmp07*transpose(G0_k); n1k=C_mat(end)*gk_7; G_07_c=n1k/(C_mat(end)*n_1k0);
gk_8=tmp08*transpose(G0_k); n1k=C_mat(end)*gk_8; G_08_c=n1k/(C_mat(end)*n_1k0);
gk_9=tmp09*transpose(G0_k); n1k=C_mat(end)*gk_9; G_09_c=n1k/(C_mat(end)*n_1k0);
gk_10=tmp10*transpose(G0_k);n1k=C_mat(end)*gk_10;G_10_c=n1k/(C_mat(end)*n_1k0);

G_f_c=[G_00_c G_01_c G_02_c G_03_c G_04_c G_05_c G_06_c G_07_c G_08_c G_09_c G_10_c];
dm_c=G_f_c./term;

%fprintf('real dm = %f \n',real(dm_c));
%fprintf('imag dm = %f \n',imag(dm_c));
%fprintf('ratio = %f \n',imag(dm_c)./real(dm_c));

fprintf('d_m = %5.5f \n',real(dm_c));

%{
figure(99); max_lim=abs(G_c(end)); x_fake=[0 max_lim 0 -max_lim]; y_fake=[max_lim 0 -max_lim 0]; h_fake=compass(x_fake,y_fake); hold on; set(h_fake,'Visible','off');
compass(G_00_c,'r'); hold on;
compass(G_01_c,'g'); hold on;
compass(G_02_c,'b'); hold on;
compass(G_03_c,'c'); hold on;
compass(G_04_c,'black'); hold on;
compass(G_05_c,'r'); hold on;
compass(G_06_c,'g'); hold on;
compass(G_c(end),'b'); hold on;
%}

pause;
