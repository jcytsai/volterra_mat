function func=dsigmadelta_ds(s_array)

global re Nb s_ele C_ele egamma_vec sigma_z0;

egamma_array=interp1(s_ele,egamma_vec,s_array);
%sigma_zs_prime=sigma_z0./interp1(s_ele,C_ele,s_array);

func=0.22*re*Nb./egamma_array;
%func=0.22*re*Nb./(egamma_array.*sigma_zs_prime.^(4/3));