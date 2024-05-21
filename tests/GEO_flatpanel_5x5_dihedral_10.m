function [AC,State] = GEO_flatpanel_5x5_dihedral_10
State.Qinf = [274.4*cosd(5) 0 274.4*sind(5)];% Velocity [Vx Vy Vz] for static loadings, m/s
State.rho_air = 1.225;%kg/m^3
State.M = 0.8;% Mach number
State.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)

%% Wing section 1 horizontal
AC(1).Label='Wing';
AC(1).Config.Proot=[0 0 0];%m  [x  y  z]
AC(1).Config.Ptip=[0 4 0];%m   [x  y  z]
AC(1).Config.rChord=1;%m
AC(1).Config.tChord=1;%m
AC(1).Config.Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).Config.Airfoil='flat';% Airfoil file
AC(1).Config.Dihedral=10*pi/180;% dihedral angle
AC(1).Config.Dihedral_ref = 1;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(1).Config.Symmetry='no';% symmetrical wing

% Aero Analysis (Gen Panel)                                                       
AC(1).SubSurf(1).Name='Wing';% wing surface
AC(1).SubSurf(1).SpanData=[0 1 5];% span percentage interval, nspan
AC(1).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(1).SubSurf(1).ChordMesh=[5];%nchord1, nchord2
AC(1).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

AcControl = [];
AcTrim = [];