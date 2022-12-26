function func=dz_dtheta(L,theta,rho,egamma)

% The input argument theta is a vector and the others are scalars.
% The input argument L can be a vector, too.

tmp01=1/(2*egamma^2);
tmp02=(theta.*(theta+2*L/rho)).^2;
tmp03=8*(theta+L/rho).^2;

func=rho*(tmp01+tmp02./tmp03);