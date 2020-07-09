function func=df_kernel_dtheta(L,theta,rho,egamma)

% The input argument theta is a vector and the others are scalars.
% The input argument L can be a vector, too.

egammatheta2=(egamma*theta).^2;
egammatheta4=(egamma*theta).^4;

tmp01=4*L.^2.*(2+egammatheta2+egammatheta4);
tmp02=4*L.*theta*rho.*(4+2*egammatheta2+egammatheta4);
tmp03=rho^2*theta.^2.*(8+2*egammatheta2+egammatheta4);
tmp04=4*L.^2.*(1+egammatheta2);
tmp05=4*L.*theta*rho.*(2+egammatheta2);
tmp06=rho^2*theta.^2.*(4+egammatheta2);

func=-rho^2*(tmp01+tmp02+tmp03)./((tmp04+tmp05+tmp06).^2);