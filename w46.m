function func=w46(s_array)

global s_ele C_ele R46_ele Nb sigma_z0;

[~,c]=size(s_array);
func=zeros(1,c);

for p=2:1:c
    s_lower=s_array(1);
    s_upper=s_array(p);
    s_prime=linspace(s_lower,s_upper,1+(p-1)*10);
    R46s_prime=interp1(s_ele,R46_ele,s_prime);
    [hBys_prime,rhoys_prime]=hBy(s_prime);
    sigma_zs_prime=sigma_z0./interp1(s_ele,C_ele,s_prime);
    X=s_prime;
    Y=R46s_prime.*hBys_prime.*sigma_zs_prime.^(-4/3).*rhoys_prime.^(-2/3);
    func(p)=trapz(X,Y);
end

