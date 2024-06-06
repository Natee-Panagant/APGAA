function [AC,FC] = GEO_further_away_case
FC.Qinf = [274.4*cosd(5) 0 274.4*sind(5)];% Velocity [Vx Vy Vz] for static loadings, m/s
FC.rho_air = 1.225;%kg/m^3
FC.M = 0.8;% Mach number
FC.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)

%% Wing section 1 horizontal
wingXpos=4.5;%wing x position, meter
wingYpos=0;%wing y position, meter
wingZpos=0;%wing z position, meter
AC(1).Label='Wing';
AC(1).Proot=[wingXpos wingYpos wingZpos];%m
AC(1).Ptip=[wingXpos wingYpos+3 wingZpos];%m
AC(1).rChord=1.5;%m
AC(1).tChord=1.5;%m
AC(1).nSpanPanel=2;% no. of spanwise vortex panels
AC(1).nChordPanel=2;% no. of chordwise vortex panels
AC(1).Dihedral=0*pi/180;% dihedral angle
AC(1).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(1).Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(1).Airfoil='flat';% Airfoil file
AC(1).Symmetry='no';% symmetrical wing

%% Winglet
wingXpos=4.5;%wing x position, meter
wingYpos=3;%wing y position, meter
wingZpos=0;%wing z position, meter
AC(2).Label='Wing';
AC(2).Proot=[wingXpos wingYpos wingZpos];%m
AC(2).Ptip=[wingXpos wingYpos wingZpos+3];%m
AC(2).rChord=1.3;%m
AC(2).tChord=1.3;%m
AC(2).nSpanPanel=2;% no. of spanwise vortex panels
AC(2).nChordPanel=2;% no. of chordwise vortex panels
AC(2).Dihedral=0*pi/180;% dihedral angle
AC(2).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(2).Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(2).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(2).Airfoil='flat';% Airfoil file
AC(2).Symmetry='no';% symmetrical wing

%%Htail
HtailXpos=11.5;%Horizontal tail x position, meter
HtailYpos=0;%Horizontal tail y position, meter
HtailZpos=1.5;%Horizontal tail z position, meter
AC(3).Label='Horizontal tail';
AC(3).Proot=[HtailXpos HtailYpos HtailZpos];%m
AC(3).Ptip=[HtailXpos HtailYpos+3 HtailZpos];%m
AC(3).rChord=1;%m
AC(3).tChord=1;%m
AC(3).nSpanPanel=2;% no. of spanwise vortex panels
AC(3).nChordPanel=2;% no. of chordwise vortex panels
AC(3).Dihedral=0*pi/180;% dihedral angle
AC(3).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(3).Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(3).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(3).Airfoil='flat';% Airfoil file
AC(3).Symmetry='no';% symmetrical wing