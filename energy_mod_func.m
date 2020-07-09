function func=energy_mod_func(s,Gs,lambda)

global egamma s_ele C_ele emitx emity betax0 betay0 sigma_delta re nb;
global iCSR_ss iCSR_tr iCSR_drift iLSC ilinac sig_x_ele sig_y_ele;

s_beg=s(1); s_end=s(end);
if (s_beg~=0); fprintf('starting not from 0...\n'); end

[r,mesh_num]=size(s); s=linspace(s_beg,s_end,mesh_num); s=s';

tmp01=re*nb/egamma;
tmp02=interp1(s_ele,C_ele,s);  % C(s')
k_wave=2*pi/lambda;
delta_p_func=zeros(mesh_num,1);

    for s_ind=1:1:mesh_num
        if (s_ind==1)
            delta_p_func(s_ind)=0.0;
        else
            tmp03=zeros(s_ind,1);
        
            for sp_ind=1:1:s_ind
                
                if ((iLSC ==0) && (ilinac ==0))
                    if ((iCSR_ss == 0)&&(iCSR_tr == 0)&&(iCSR_drift == 0))
                        tmp03(sp_ind,1)=0.0;
                    elseif ((iCSR_ss == 1)&&(iCSR_tr == 0)&&(iCSR_drift == 0))
                        tmp03(sp_ind,1)=csr1d(k_wave*tmp02(sp_ind),s(sp_ind));
                    elseif ((iCSR_ss == 0)&&(iCSR_tr == 1)&&(iCSR_drift == 0))
                        tmp03(sp_ind,1)=csr1d_tr(k_wave*tmp02(sp_ind),s(sp_ind));
                    elseif ((iCSR_ss == 0)&&(iCSR_tr == 0)&&(iCSR_drift == 1))
                        tmp03(sp_ind,1)=csr1d_drift(k_wave*tmp02(sp_ind),s(sp_ind));
                    elseif ((iCSR_ss == 1)&&(iCSR_tr == 1)&&(iCSR_drift == 0))
                        tmp03(sp_ind,1)=csr1d(k_wave*tmp02(sp_ind),s(sp_ind))+csr1d_tr(k_wave*tmp02(sp_ind),s(sp_ind));
                    elseif ((iCSR_ss == 1)&&(iCSR_tr == 0)&&(iCSR_drift == 1))
                        tmp03(sp_ind,1)=csr1d(k_wave*tmp02(sp_ind),s(sp_ind))+csr1d_drift(k_wave*tmp02(sp_ind),s(sp_ind));
                    elseif ((iCSR_ss == 0)&&(iCSR_tr == 1)&&(iCSR_drift == 1))
                        tmp03(sp_ind,1)=csr1d_tr(k_wave*tmp02(sp_ind),s(sp_ind))+csr1d_drift(k_wave*tmp02(sp_ind),s(sp_ind));
                    elseif ((iCSR_ss == 1)&&(iCSR_tr == 1)&&(iCSR_drift == 1))
                        tmp03(sp_ind,1)=csr1d(k_wave*tmp02(sp_ind),s(sp_ind))+csr1d_tr(k_wave*tmp02(sp_ind),s(sp_ind))+csr1d_drift(k_wave*tmp02(sp_ind),s(sp_ind));
                    end
                end
                
                % LSC-related models
                if (iLSC==1)
                    tmp_LSC=lsc1d(k_wave*tmp02(sp_ind),interp1(s_ele,sig_x_ele,s(sp_ind)),interp1(s_ele,sig_y_ele,s(sp_ind)));
                else
                    tmp_LSC=0.0;
                end
    
                tmp03(sp_ind,1)=tmp03(sp_ind,1)+tmp_LSC;

            end
            sp=s(1:s_ind);
            
            tmp06=0.5*k_wave^2*emitx/(betax0);
            tmp07=betax0^2*R51mod(s(s_ind),sp).^2+R52mod(s(s_ind),sp).^2; % vector
            tmp08=0.5*(k_wave^2)*(sigma_delta^2)*(R56mod(s(s_ind),sp).^2); % vector
            tmp09=0.5*k_wave^2*emity/(betay0);
            tmp10=betay0^2*R53mod(s(s_ind),sp).^2+R54mod(s(s_ind),sp).^2; % vector
            tmp04=Gs(1:s_ind);
            tmp05=exp(-tmp06.*tmp07-tmp09.*tmp10-tmp08);
            delta_p_func(s_ind)=abs(trapz(s(1:s_ind),-tmp01*tmp02(1:s_ind).*tmp03.*tmp04.*tmp05));
        end
    end

func=delta_p_func;



