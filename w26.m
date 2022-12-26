function func=w26(s_array)

global s_ele C_ele R26_ele Nb sigma_z0;

[~,c]=size(s_array);
func=zeros(1,c);

for p=2:1:c
    s_lower=s_array(1);
    s_upper=s_array(p);
    s_prime=linspace(s_lower,s_upper,1+(p-1)*10);
    R26s_prime=interp1(s_ele,R26_ele,s_prime);
    [hBxs_prime,rhoxs_prime]=hBx(s_prime);
    sigma_zs_prime=sigma_z0./interp1(s_ele,C_ele,s_prime);
    X=s_prime;
    Y=R26s_prime.*hBxs_prime.*sigma_zs_prime.^(-4/3).*rhoxs_prime.^(-2/3);
    %Y=R26s_prime.*hBxs_prime.*rhoxs_prime.^(-2/3);
    func(p)=trapz(X,Y);
end

