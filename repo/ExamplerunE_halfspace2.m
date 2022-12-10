%Example to run for the site results
%https://em.geosci.xyz/content/maxwell3_fdem/inductive_sources/halfspace.html
%Figure 104,(a)
clear all
% Transmitter coordinates
xT=00;
yT=00;
zT=-30;
%Conductivity
cond=0.01;
%Frequency
freq=10e3;


% For specific y
yy=0;
countx=0;
xvector=-100.001:1:99.999; % force it not to go to zero
for xx=xvector
    countx=countx+1;
    countz=0;
    zvector=-100:1:100;
        for zz=zvector
            countz=countz+1;
            AB=[xx,yy,zz];

[Ey, Ex]=E_to_onepoint_hankel_transform(AB, xT, yT, zT, cond, freq);

Eygrid(countz,countx)=Ey;

Exgrid(countz,countx)=Ex;
        end
end

% Plot
figure;imagesc(xvector,zvector,real(Eygrid));colorbar

xlabel('Distance (m)'); ylabel('Depth (m)');
 title(['Ey calculated, zT=',num2str(zT),'m, sigma=', num2str(cond),'S/m, Frequency=', num2str(freq),'Hz'])
 hold on
 plot(xT,zT,'ob')
