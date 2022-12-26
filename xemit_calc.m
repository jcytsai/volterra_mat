function xemit_calc(s_array)

global s_ele_Twiss alphax_s_ele alphay_s_ele betax_s_ele betay_s_ele gammax_s_ele gammay_s_ele Nb;
global s_ele R16_ele R26_ele R36_ele R46_ele egamma_vec emitx emity re;

tmp01=Dxp2(s_array);
tmp02=DxDxp(s_array);
tmp03=Dx2(s_array);
tmp04=tmp03.*tmp01;
tmp05=tmp02.^2;

alphax=interp1(s_ele_Twiss,alphax_s_ele,s_array);
betax=interp1(s_ele_Twiss,betax_s_ele,s_array);
gammax=interp1(s_ele_Twiss,gammax_s_ele,s_array);

emitx_CSR=sqrt(emitx^2+emitx*(betax.*tmp01+2*alphax.*tmp02+gammax.*tmp03)+tmp04-tmp05);

figure(401); set(gca,'FontSize',40,'linewidth',5); plot(s_array/100,emitx_CSR/emitx,'r-','linewidth',5);   xlabel('s (m)'); ylabel('\epsilon_x^{CSR}/\epsilon_{x0}'); grid off; hold on; axis('tight');

