function[Ex,Ey]=magnetic_dipole_N(f,source,receiver,sigma)
omega=2*pi*f; mi=4*pi*10^-7; % f =30000; sigma=10^-8 or 10^-2;
m=1;k=sqrt(-1i*mi*sigma*omega); %m: magnetic strength 
x=receiver(1,:)-source(1,:);
y=receiver(2,:)-source(2,:);
z=receiver(3,:)-source(3,:);
r=(x.*x+y.*y+z.*z).^0.5; %T-R distance
Efactor=1i*omega*mi*m/4/pi./r.^3.*(1i*k.*r +1).*exp(-1i*k.*r); % Eqtn 2.56 for M vertical front
E=[Efactor.*y; -Efactor.*x; 0.*z]; % Eqtn 2. 56 modified E for z-direction magnetic dipole in whole space
Ex=E(1,:);
Ey=E(2,:);
end