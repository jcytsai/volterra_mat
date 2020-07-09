function func=energy_mod_spec(s,G_matrix,lambda_start,scan_mesh)

global egamma s_ele C_ele emitx emity betax0 betay0 sigma_delta re nb;
global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac sig_x_ele sig_y_ele;

s_beg=s(1); s_end=s(end);
if (s_beg~=0); fprintf('starting not from 0...\n'); end

[r,mesh_num]=size(s); s=linspace(s_beg,s_end,mesh_num); s=s';
[r,scan_num]=size(G_matrix); % row: s, col: iscan

tmp01=re*nb/egamma;
tmp02=interp1(s_ele,C_ele,s);   % C(s')

delta_p_spec=zeros(1,scan_num);
for iscan=1:1:scan_num
    lambda=lambda_start+(iscan-1)*scan_mesh;
    k_wave=2*pi/lambda;
    
    tmp03=zeros(mesh_num,1);
    tmp06=zeros(mesh_num,1);
    tmp07=zeros(mesh_num,1);
    tmp08=zeros(mesh_num,1);
    tmp09=zeros(mesh_num,1);
    tmp10=zeros(mesh_num,1);
    
    for s_ind=1:1:mesh_num
        sp=s(s_ind,1);
        
        if ((iLSC ==0) && (ilinac ==0))
            if ((iCSR_ss == 0)&&(iCSR_tr == 0)&&(iCSR_drift == 0))
                tmp03(s_ind,1)=0.0;
            elseif ((iCSR_ss == 1)&&(iCSR_tr == 0)&&(iCSR_drift == 0))
                tmp03(s_ind,1)=csr1d(k_wave*tmp02(s_ind,1),sp);
            elseif ((iCSR_ss == 0)&&(iCSR_tr == 1)&&(iCSR_drift == 0))
                tmp03(s_ind,1)=csr1d_tr(k_wave*tmp02(s_ind,1),sp);
            elseif ((iCSR_ss == 0)&&(iCSR_tr == 0)&&(iCSR_drift == 1))
                tmp03(s_ind,1)=csr1d_drift(k_wave*tmp02(s_ind,1),sp);
            elseif ((iCSR_ss == 1)&&(iCSR_tr == 1)&&(iCSR_drift == 0))
                tmp03(s_ind,1)=csr1d(k_wave*tmp02(s_ind,1),sp)+csr1d_tr(k_wave*tmp02(s_ind,1),sp);
            elseif ((iCSR_ss == 1)&&(iCSR_tr == 0)&&(iCSR_drift == 1))
                tmp03(s_ind,1)=csr1d(k_wave*tmp02(s_ind,1),sp)+csr1d_drift(k_wave*tmp02(s_ind,1),sp);
            elseif ((iCSR_ss == 0)&&(iCSR_tr == 1)&&(iCSR_drift == 1))
                tmp03(s_ind,1)=csr1d_tr(k_wave*tmp02(s_ind,1),sp)+csr1d_drift(k_wave*tmp02(s_ind,1),sp);
            elseif ((iCSR_ss == 1)&&(iCSR_tr == 1)&&(iCSR_drift == 1))
                tmp03(s_ind,1)=csr1d(k_wave*tmp02(s_ind,1),sp)+csr1d_tr(k_wave*tmp02(s_ind,1),sp)+csr1d_drift(k_wave*tmp02(s_ind,1),sp);
            end
        end
        
        % LSC-related models
        if (iLSC==1)
        	tmp_LSC=lsc1d(k_wave*tmp02(s_ind,1),interp1(s_ele,sig_x_ele,sp),interp1(s_ele,sig_y_ele,sp));
        else
            tmp_LSC=0.0;
        end
        tmp03(sp_ind,1)=tmp03(sp_ind,1)+tmp_LSC;
        
        tmp06(s_ind,1)=0.5*k_wave^2*emitx/(betax0);
        tmp07(s_ind,1)=betax0^2*R51mod(s(end,1),sp).^2+R52mod(s(end,1),sp).^2; % vector
        tmp08(s_ind,1)=0.5*(k_wave^2)*(sigma_delta^2)*(R56mod(s(end,1),sp).^2); % vector
        tmp09(s_ind,1)=0.5*k_wave^2*emity/(betay0);
        tmp10(s_ind,1)=betay0^2*R53mod(s(end,1),sp).^2+R54mod(s(end,1),sp).^2; % vector       
    end
    
    tmp04=G_matrix(:,iscan);
    tmp05=exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08);
    
    delta_p_spec(1,iscan)=abs(trapz(s,-tmp01*tmp02.*tmp03.*tmp04.*tmp05));
    
end

func=delta_p_spec;




