%clear all
format long

c_speed=2.99792458e8; % m/s

phi_deg_select=[-12 -11 -10 -9 -8 8 9 10 11 12]; % deg
phi_select=phi_deg_select/180*pi;

phi_deg=linspace(-15,15,3000); % deg
phi=phi_deg/180*pi;

E0=1.3e3; % MeV
frf=1497e6; % Hz
k=2*pi/(c_speed/frf);

eV=1.2e3; % MeV

R65=eV/E0*k*sin(phi);
T655=-0.5*eV/E0*k^2*cos(phi);

R65_select=eV/E0*k*sin(phi_select);
T655_select=-0.5*eV/E0*k^2*cos(phi_select);

figure(1); plot(phi_deg,R65); xlabel('\phi (deg)'); ylabel('R_{65} (m^{-1})'); hold on;
figure(1); plot(phi_deg_select,R65_select,'o'); xlabel('\phi (deg)'); ylabel('R_{65} (m^{-1})'); hold on;

figure(2); plot(phi_deg,T655); xlabel('\phi (deg)'); ylabel('T_{655} (m^{-2})'); hold on;
figure(2); plot(phi_deg_select,T655_select,'o'); xlabel('\phi (deg)'); ylabel('T_{655} (m^{-2})'); hold on;