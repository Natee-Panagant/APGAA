function [AC,FC] = GEO_flatpanel_5x5_dihedral_10
FC.Qinf = [274.4*cosd(5) 0 274.4*sind(5)];% Velocity [Vx Vy Vz] for static loadings, m/s
FC.rho_air = 1.225;%kg/m^3
FC.M = 0.8;% Mach number
FC.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)
FC.rG = [0 0 0];%Reference point for moment caculation

%% Wing section 1 horizontal
wingXpos=0;%wing x position, meter
wingYpos=0;%wing y position, meter
wingZpos=0;%wing z position, meter
AC(1).Label='Wing';
AC(1).Proot=[wingXpos wingYpos wingZpos];%m
AC(1).Ptip=[wingXpos wingYpos+4 wingZpos];%m
AC(1).rChord=1;%m
AC(1).tChord=1;%m
AC(1).nSpanPanel=5;% no. of spanwise vortex panels
AC(1).nChordPanel=5;% no. of chordwise vortex panels
AC(1).Dihedral=10*pi/180;% dihedral angle
AC(1).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(1).Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(1).Airfoil='flat';% Airfoil file
AC(1).Symmetry='no';% symmetrical wing