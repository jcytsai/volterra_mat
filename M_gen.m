function [func1,func2]=M_gen(s_shift,section_length,egamma0,RF_freq,psi0_deg,E0,M_ini)

global alpha_RF lambda_RF k_RF omega_RF psi0;

c_speed=2.99792458e8;                  % unit: m/sec

t0=s_shift/c_speed;
s0=s_shift;

M11_ini=M_ini(1,1);
M12_ini=M_ini(1,2);
M13_ini=M_ini(1,3);
M14_ini=M_ini(1,4);
M15_ini=M_ini(1,5);
M16_ini=M_ini(1,6);

M21_ini=M_ini(2,1);
M22_ini=M_ini(2,2);
M23_ini=M_ini(2,3);
M24_ini=M_ini(2,4);
M25_ini=M_ini(2,5);
M26_ini=M_ini(2,6);

M31_ini=M_ini(3,1);
M32_ini=M_ini(3,2);
M33_ini=M_ini(3,3);
M34_ini=M_ini(3,4);
M35_ini=M_ini(3,5);
M36_ini=M_ini(3,6);

M41_ini=M_ini(4,1);
M42_ini=M_ini(4,2);
M43_ini=M_ini(4,3);
M44_ini=M_ini(4,4);
M45_ini=M_ini(4,5);
M46_ini=M_ini(4,6);

M51_ini=M_ini(5,1);
M52_ini=M_ini(5,2);
M53_ini=M_ini(5,3);
M54_ini=M_ini(5,4);
M55_ini=M_ini(5,5);
M56_ini=M_ini(5,6);

M61_ini=M_ini(6,1);
M62_ini=M_ini(6,2);
M63_ini=M_ini(6,3);
M64_ini=M_ini(6,4);
M65_ini=M_ini(6,5);
M66_ini=M_ini(6,6);


psi0=psi0_deg/180*pi;                  % unit: rad
lambda_RF=c_speed/RF_freq;             % unit: m
beta_rel=sqrt(1-1/egamma0^2);
omega_RF=2*pi*RF_freq;
k_RF=2*pi/beta_rel/lambda_RF;
alpha_RF=E0/(0.511*k_RF);

Y0=[t0;egamma0;...
    M11_ini;M12_ini;M13_ini;M14_ini;M15_ini;M16_ini;...
    M21_ini;M22_ini;M23_ini;M24_ini;M25_ini;M26_ini;...
    M31_ini;M32_ini;M33_ini;M34_ini;M35_ini;M36_ini;...
    M41_ini;M42_ini;M43_ini;M44_ini;M45_ini;M46_ini;...
    M51_ini;M52_ini;M53_ini;M54_ini;M55_ini;M56_ini;...
    M61_ini;M62_ini;M63_ini;M64_ini;M65_ini;M66_ini;...
    s0];

s_span=[s0 s0+section_length];

[~,Y]=ode45(@RHS_eq,s_span,Y0,'options',E0); 
options=odeset('RelTol',1e-10,'AbsTol',1e-10,'Stats','on');

t=Y(:,1);                              % time of flight
egamma_vec=Y(:,2);                     % egamma
R11_vec=Y(:,3);                        % R11(s)
R12_vec=Y(:,4);                        % R12(s)
R13_vec=Y(:,5);                        % R13(s)
R14_vec=Y(:,6);                        % R14(s)
R15_vec=Y(:,7);                        % R15(s)
R16_vec=Y(:,8);                        % R16(s)
R21_vec=Y(:,9);                        % R21(s)
R22_vec=Y(:,10);                       % R22(s)
R23_vec=Y(:,11);                       % R23(s)
R24_vec=Y(:,12);                       % R24(s)
R25_vec=Y(:,13);                       % R25(s)
R26_vec=Y(:,14);                       % R26(s)
R31_vec=Y(:,15);                       % R31(s)
R32_vec=Y(:,16);                       % R32(s)
R33_vec=Y(:,17);                       % R33(s)
R34_vec=Y(:,18);                       % R34(s)
R35_vec=Y(:,19);                       % R35(s)
R36_vec=Y(:,20);                       % R36(s)
R41_vec=Y(:,21);                       % R41(s)
R42_vec=Y(:,22);                       % R42(s)
R43_vec=Y(:,23);                       % R43(s)
R44_vec=Y(:,24);                       % R44(s)
R45_vec=Y(:,25);                       % R45(s)
R46_vec=Y(:,26);                       % R46(s)
R51_vec=Y(:,27);                       % R51(s)
R52_vec=Y(:,28);                       % R52(s)
R53_vec=Y(:,29);                       % R53(s)
R54_vec=Y(:,30);                       % R54(s)
R55_vec=Y(:,31);                       % R55(s)
R56_vec=Y(:,32);                       % R56(s)
R61_vec=Y(:,33);                       % R61(s)
R62_vec=Y(:,34);                       % R62(s)
R63_vec=Y(:,35);                       % R63(s)
R64_vec=Y(:,36);                       % R64(s)
R65_vec=Y(:,37);                       % R65(s)
R66_vec=Y(:,38);                       % R66(s)
s=Y(:,39);                             % s

egamma=egamma_vec(end);

func1=egamma_vec;
func2=[s R11_vec R12_vec R13_vec R14_vec R15_vec R16_vec ...
         R21_vec R22_vec R23_vec R24_vec R25_vec R26_vec ...
         R31_vec R32_vec R33_vec R34_vec R35_vec R36_vec ...
         R41_vec R42_vec R43_vec R44_vec R45_vec R46_vec ...
         R51_vec R52_vec R53_vec R54_vec R55_vec R56_vec ...
         R61_vec R62_vec R63_vec R64_vec R65_vec R66_vec];
         
         
         



