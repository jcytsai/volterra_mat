function func=dz(L,theta,rho,egamma)

% The input argument theta is a vector and the others are scalars.
% The input argument L can be a vector, too.

ebeta=sqrt(1-egamma^(-2));
tmp01=(L+rho*theta)/(2*egamma^2);
tmp02=ebeta*rho*theta.^3/24;
tmp03=(4*L+rho*theta)./(L+rho*theta);

func=tmp01+tmp02.*tmp03;