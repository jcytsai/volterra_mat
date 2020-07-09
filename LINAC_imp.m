function func=LINAC_imp(k,a,L,g,alpha,s)

global RF_info;

% k: scaled wave number in unit of 1/cm
% a: average iris radius in cm
% L: cell(period) length in cm
% g: gap distance between irises
% alpha: function of g/L (~0.5 for g/L~1)
% s: position in cm, used to determine if the impedance is to be added

%{
tmp01=find((RF_info(:,1)-s/100)<0,1,'last');

if (mod(tmp01,2)==0)
    tmp04=0.0;
else
    tmp02=4*1i/(k*a^2);
    tmp03=1+(1+1i)*alpha*L/a*sqrt(pi/g/k);
    tmp04=tmp02/tmp03;
    %fprintf('linac geometric impedance Z=%f+%fi applied at s = %f m...\n',real(tmp04),imag(tmp04),s/100);
end

func=tmp04;
%}

tmp02=4*1i/(k*a^2);
tmp03=1+(1+1i)*alpha*L/a*sqrt(pi/g/k);
tmp04=tmp02/tmp03;
%fprintf('linac geometric impedance Z=%f+%fi applied at s = %f m...\n',real(tmp04),imag(tmp04),s/100);

func=tmp04;

