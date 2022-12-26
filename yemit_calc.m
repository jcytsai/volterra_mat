function yemit_calc(s_array)

global s_ele_Twiss alphax_s_ele alphay_s_ele betax_s_ele betay_s_ele gammax_s_ele gammay_s_ele Nb;
global s_ele R16_ele R26_ele R36_ele R46_ele egamma_vec emitx emity re;

tmp01=Dyp2(s_array);
tmp02=DyDyp(s_array);
tmp03=Dy2(s_array);
tmp04=tmp03.*tmp01;
tmp05=tmp02.^2;

alphay=interp1(s_ele_Twiss,alphay_s_ele,s_array);
betay=interp1(s_ele_Twiss,betay_s_ele,s_array);
gammay=interp1(s_ele_Twiss,gammay_s_ele,s_array);

emity_CSR=sqrt(emity^2+emity*(betay.*tmp01+2*alphay.*tmp02+gammay.*tmp03)+tmp04-tmp05);

figure(401); set(gca,'FontSize',40,'linewidth',5); plot(s_array/100,emity_CSR/emity,'b-','linewidth',5);   xlabel('s (m)'); ylabel('\epsilon_y^{CSR}/\epsilon_{y0}'); grid off; hold on; axis('tight');

