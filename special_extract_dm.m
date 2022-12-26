% this special program intends to extract stage gain coefficient d_m
% directly from G_f^(M), instead of using polynominal fitting technique.
% For more information, see Eqs.(10,11,14) of the long paper.
% To run this subroutine, one needs to run a case with lambda_opt.

Alfven=17045;

max_order=10;
I_b_start=480;
I_b_end=480;
I_b_num=1;
I_b_array=linspace(I_b_start,I_b_end,I_b_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coeff=const_A*I_b/(egamma*Alfven);
coeff_array=const_A*I_b_array/(egamma*Alfven);
term=coeff.^(0:1:max_order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp01=mesh_size*K_mat(end,:);   
tmp02=tmp01*(mesh_size*K_mat); 
tmp03=tmp02*(mesh_size*K_mat); 
tmp04=tmp03*(mesh_size*K_mat); 
tmp05=tmp04*(mesh_size*K_mat); 
tmp06=tmp05*(mesh_size*K_mat); 
tmp07=tmp06*(mesh_size*K_mat); 
tmp08=tmp07*(mesh_size*K_mat); 
tmp09=tmp08*(mesh_size*K_mat); 
tmp10=tmp09*(mesh_size*K_mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
C_mat=interp1(s_ele,C_ele,s); C_mat=transpose(C_mat);

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
%
IAmdm_c=dm_c./((Alfven).^(0:1:max_order));
%fprintf('real IAmdm = %f \n',real(IAmdm_c));
%fprintf('imag IAmdm = %f \n',imag(IAmdm_c));
fprintf('abs IAmdm = %f \n',abs(IAmdm_c));
%}
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

stage_gain_c_00=dm_c(1)*coeff_array.^0;
stage_gain_c_01=dm_c(2)*coeff_array.^1;
stage_gain_c_02=dm_c(3)*coeff_array.^2;
stage_gain_c_03=dm_c(4)*coeff_array.^3;
stage_gain_c_04=dm_c(5)*coeff_array.^4;
stage_gain_c_05=dm_c(6)*coeff_array.^5;
stage_gain_c_06=dm_c(7)*coeff_array.^6;
stage_gain_c_07=dm_c(8)*coeff_array.^7;
stage_gain_c_08=dm_c(9)*coeff_array.^8;
stage_gain_c_09=dm_c(10)*coeff_array.^9;
stage_gain_c_10=dm_c(11)*coeff_array.^10;

%
stage_gain_c=stage_gain_c_00+stage_gain_c_01+stage_gain_c_02+...
             stage_gain_c_03+stage_gain_c_04+stage_gain_c_05+...
             stage_gain_c_06+stage_gain_c_07+stage_gain_c_08+...
             stage_gain_c_09+stage_gain_c_10;
%}
%{
stage_gain_c=stage_gain_c_00+stage_gain_c_01+stage_gain_c_02+...
             stage_gain_c_03+stage_gain_c_04+stage_gain_c_05+...
             stage_gain_c_06;
%}
figure(999); set(gca,'FontSize',40,'linewidth',5); semilogy(I_b_array,abs(stage_gain_c)); xlabel('I_b (A)'); ylabel('gain'); hold on;




