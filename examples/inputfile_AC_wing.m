function [AC,FC] = inputfile_AC_wing

%% 

%%
% Flight conditions
FC.Qinf = [10 0 0];% Velocity for static loadings, m/s
FC.rho_air = 1.20;% kg/m^3
FC.M = 0.5;% Mach number
FC.rG = [1.00 0 0];% center of mass
FC.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)


%%
% Aircraft configuration details

%%
% Wing
wingXpos = 3; %wing x position, meter
wingYpos = 0.6;%wing y position, meter
wingZpos = -0.4;%wing z position, meter

AC(1).Label = 'Wing';
AC(1).Proot = [wingXpos wingYpos wingZpos];%m
AC(1).Ptip = [wingXpos+0.4 wingYpos+10 wingZpos];%m
AC(1).rChord = 1.25;    %m
AC(1).tChord = 1.2;      %m
AC(1).nSpanPanel = 20;% no. of spanwise vortex panels
AC(1).nChordPanel = 5;% no. of chordwise vortex panels
AC(1).Dihedral = 5*pi/180;% dihedral angle
AC(1).dRefPt = [0 0 wingZpos];% % reference point for dihedral angle rotating
AC(1).Incidence = 3*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).iRefPt = [wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(1).Airfoil = 'flat';% Airfoil file
AC(1).Symmetry = 'yes';% symmetrical wing  


end
