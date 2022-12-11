function [Hhank, Hn] = inifile_econ1(AB, xT, yT, zT, cond, freq, orient) 

% Info
% Calculates electric or magnetic fields either in the Earth (only the electric) or above it,
% with a transmitter loop placed either on the Earth's surface of above it
% (this is why zT<=0 always).

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- User-input -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AB is a vector containing [xR, yR, zR], which are the coordinates of the transmitter in meters. zR is positive below the Earth's surface and negative above it.
% xT is the x-coordinate of the transmitter (m)
% yT is the y-coordinate of the transmitter (m)
% zT is the z-coordinate of the transmitter (m)-Use a negative number only,which means you always have the transmitter above the Earth
% cond is the conductivity in S/m
% freq is the frequency in Herz
% orient is a text label indicating whether you want the magnetic field of a vertical magnetic dipole ('ZZ'), or the electric field of a vertical magnetic dipole ('ZZE')

% -------------------------------------------------------Output-----------%
% Hhank is the total field response. Its Ef in cylidrical coordinates and must be transformed in cartesian coordinates, if orient is 'ZZE'. If
% it is 'ZZ' then we have Hz directly.


% Hn is the air field response -primary-background (Hzhank-Hn=Hres, which
% is the secondary fied). The same as above for the coordinates.


% If orient = 'ZZ', then we have the magnetic fields, if orient is 'ZZE' then we have the electric fields.

% If you choose the estimation of the magnetic field dont bother estimating
% below the Earth's surface beacuse this program sets it to zero if the
% receiver is below the Earth (line 168 in FDEM1DFWD_RC1). For the electric fields you can choose a
% receiver position wherever you want (try to avoid the Transmitter's
% position for both cases of electric or magnetic fields).

%------------------------------------------------------Example------------%

%[Hhank, Hn] = inifile_econ1([1 0 -1], 0, 0, -1, 0.01, 10e3, 'ZZE') 
% The above calculates the Electric fields at a point at x-coord =1,
% y-coord=0 and 1 meter above the Earth's surface, while the Transmitter is
% at x-coord=0, y-coord=0 and is placed 1 meter above the Earths's surface.

%--------------------------------------------------------------------------
%%

xR=AB(1);yR=AB(2);zR=AB(3);
S.x=xR-xT;                            % x-coordinate receiver (m)
S.y=yR-yT;                            % y-coordinate receiver (m)
S.z=zR;                               % z-coordinate receiver (m) - positive z-axis pointed down
S.height=-zT;
S.mom=1;
S.freq=freq;
S.ori=orient;

    %
    % Model characteristics (M structure) 
    % *********MADE FOR NON MAGNETIC HALF SPACE WITH CONDUCTIVITY=cond S/m.***********
    %

        M.sus=      linspace(4*pi,4*pi,1) .* 1e-7;                                                                 % Susceptibility of layer(s) (-)
        M.con=      linspace(cond,cond,1);                                                                         % Conductivity of layer(s) (S/m)
        M.perm=     linspace(1/(35950207149.4727056*pi),1/(35950207149.4727056*pi),1);                             % Permittivity of layer(s) (F/m)
        M.thick=    logspace(log10(0.1),log10(0.1),1);                                                             % Layer(s) thickness (m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- Modelling --------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % Calculate forward response (ppm)
    %

        [FWD_IP,FWD_QP,Hn] = FDEM1DFWD_RC1(S,M);                               % Reflection coefficient approach
        
Hhank=FWD_IP+i*FWD_QP;


