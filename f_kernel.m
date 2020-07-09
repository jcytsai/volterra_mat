function func=f_kernel(L,theta,rho,egamma)

% The input argument theta is a vector and the others are scalars.
% The input argument L can be a vector, too.

tmp01=2*(L/rho+theta)/egamma^2;
tmp02=(theta.^2).*(theta+2*L/rho);
tmp03=(2*(L/rho+theta)/egamma).^2;
tmp04=(theta.*(theta+2*L/rho)).^2;

func=(tmp01+tmp02)./(tmp03+tmp04);
%func=1./(theta+2*L/rho);                % ultra-relativistic approximation