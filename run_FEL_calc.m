
global MX und_lambda und_L lambda_r P_L bk_ini Pk_ini
global I_b_fin sigma_delta_fin sig_x emit_nx_fin und_beta_ave

mc2=0.511;                                          % unit: MeV

MX=0;
und_lambda=0.035;                                   % unit: m
und_L=6*68*und_lambda;                              % unit: m
lambda_r=7.3e-9;                                    % unit: m
P_L=10;                                             % seed power, unit: MW
bk_ini=3e-4;                                        % ini. shot-noise dmod amp.
Pk_ini=1e-5;                                        % ini. shot-noise emod amp.

I_b_fin=I_b*C_ele(end);                             % unit: Amp
sigma_delta_fin=SES_fin/(egamma_vec(end)*mc2);
sig_x=sig_x_ele(end)/100;                           % unit: m
emit_nx_fin=emitx/100*egamma;                       % unit: m
if (iIBS==1)
    emit_nx_fin=emit_gx_IBS/100*egamma;
end
und_beta_ave=sig_x^2/(emit_nx_fin/egamma_vec(end));

GUI_FEL_calc;