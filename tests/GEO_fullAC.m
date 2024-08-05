function [AC,FC] = GEO_fullAC
FC.Qinf = [274.4*cosd(5) 0 274.4*sind(5)];% Velocity [Vx Vy Vz] for static loadings, m/s
FC.rho_air = 1.225;%kg/m^3
FC.M = 0.8;% Mach number
FC.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)
FC.rG = [0 0 0];%Reference point for moment caculation

%% Wing section 1 horizontal
wingXpos=21;%wing x position, meter
wingYpos=2;%wing y position, meter
wingZpos=0;%wing z position, meter
AC(1).Label='Wing';
AC(1).Proot=[wingXpos wingYpos wingZpos];%m
AC(1).Ptip=[wingXpos-4 wingYpos+13 wingZpos];%m
AC(1).rChord=5;%m
AC(1).tChord=2;%m
AC(1).nSpanPanel=10;% no. of spanwise vortex panels
AC(1).nChordPanel=4;% no. of chordwise vortex panels
AC(1).Dihedral=0*pi/180;% dihedral angle
AC(1).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(1).Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(1).Airfoil='flat';% Airfoil file
AC(1).Symmetry='yes';% symmetrical wing


%% Vertical tail
VtailXpos=32;%Vertical tail x position, meter
VtailYpos=0;%Vertical tail x position, meter
VtailZpos=2;%Vertical tail z position, meter
AC(2).Label='Vertical tail';
AC(2).Proot=[VtailXpos VtailYpos VtailZpos];%m
AC(2).Ptip=[VtailXpos+5 VtailYpos VtailZpos+5];%m
AC(2).rChord=6;%m
AC(2).tChord=3;%m
AC(2).nSpanPanel=8;% no. of spanwise vortex panels
AC(2).nChordPanel=4;% no. of chordwise vortex panels
AC(2).Dihedral=0*pi/180;% dihedral angle
AC(2).dRefPt=[VtailXpos VtailYpos VtailZpos];% % reference point for dihedral angle rotating
AC(2).Incidence=0*pi/180;% incidence, radian
AC(2).iRefPt=[VtailXpos VtailYpos VtailZpos];% reference point for incidence angle rotating
AC(2).Airfoil='flat';% Airfoil file
AC(2).Symmetry='no';% symmetrical wing


%% Horizontal tail
HtailXpos=37;%Horizontal tail x position, meter
HtailYpos=0;%Horizontal tail y position, meter
HtailZpos=7;%Horizontal tail z position, meter
AC(3).Label='Horizontal tail';
AC(3).Proot=[HtailXpos HtailYpos HtailZpos];%m
AC(3).Ptip=[HtailXpos+3 HtailYpos+7 HtailZpos];%m
AC(3).rChord=3;%m
AC(3).tChord=2;%m
AC(3).nSpanPanel=8;% no. of spanwise vortex panels
AC(3).nChordPanel=4;% no. of chordwise vortex panels
AC(3).Dihedral=0*pi/180;% dihedral angle
AC(3).dRefPt=[0 0 HtailZpos];% % reference point for dihedral angle rotating
AC(3).Incidence=0*pi/180;% incidence, radian
AC(3).iRefPt=[HtailXpos HtailYpos HtailZpos];% reference point for incidence angle rotating
AC(3).Airfoil='flat';% Airfoil file
AC(3).Symmetry='yes';% symmetrical wing

%% Fuselage
FuselageXpos=0;%Fuselage tail x position, meter
FuselageYpos=0;%Fuselage tail y position, meter
FuselageZpos=0;%Fuselage z position, meter
AC(4).Label='Fuselage';
AC(4).Proot=[FuselageXpos FuselageYpos FuselageZpos];%m
AC(4).Ptip=[FuselageXpos FuselageYpos+2 FuselageZpos];%m
AC(4).rChord=40;%m
AC(4).tChord=40;%m
AC(4).nSpanPanel=1;% no. of spanwise vortex panels
AC(4).nChordPanel=20;% no. of chordwise vortex panels
AC(4).Dihedral=0*pi/180;% dihedral angle
AC(4).dRefPt=[0 0 FuselageZpos];% % reference point for dihedral angle rotating
AC(4).Incidence=0*pi/180;% incidence, radian
AC(4).iRefPt=[FuselageXpos FuselageYpos FuselageZpos];% % reference point for dihedral angle rotating
AC(4).Airfoil='flat';% Airfoil file
AC(4).Symmetry='yes';% symmetrical wing