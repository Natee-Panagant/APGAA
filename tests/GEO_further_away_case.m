function [AC,State] = GEO_further_away_case
State.Qinf = [274.4*cosd(5) 0 274.4*sind(5)];% Velocity [Vx Vy Vz] for static loadings, m/s
State.rho_air = 1.225;%kg/m^3
State.M = 0.8;% Mach number
State.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)

%% Wing section 1 horizontal
AC(1).Label='Wing';
AC(1).Config.Proot=[4.5 0 0];%m  [x  y  z]
AC(1).Config.Ptip=[4.5 3 0];%m   [x  y  z]
AC(1).Config.rChord=1.5;%m
AC(1).Config.tChord=1.5;%m
AC(1).Config.Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).Config.Airfoil='flat';% Airfoil file
AC(1).Config.Dihedral=0*pi/180;% dihedral angle
AC(1).Config.Dihedral_ref = 2;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(1).Config.Symmetry='no';% symmetrical wing
                                                   
AC(1).SubSurf(1).Name='Wing';% wing surface
AC(1).SubSurf(1).SpanData=[0 1 2];% span percentage interval, nspan
AC(1).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(1).SubSurf(1).ChordMesh=[2];%nchord1, nchord2
AC(1).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%% Winglet
AC(2).Config.Proot=[4.5 3 0];%m  [x  y  z]
AC(2).Config.Ptip=[4.5 3 3];%m   [x  y  z]
AC(2).Config.rChord=1.3;%m
AC(2).Config.tChord=1.3;%m
AC(2).Config.Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(2).Config.Airfoil='flat';% Airfoil file
AC(2).Config.Dihedral=0*pi/180;% dihedral angle
AC(2).Config.Dihedral_ref = 2;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(2).Config.Symmetry='no';% symmetrical wing

% Aero Analysis (Gen Panel)                                                       
AC(2).SubSurf(1).Name='Wing';% wing surface
AC(2).SubSurf(1).SpanData=[0 1 2];% span percentage interval, nspan
AC(2).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(2).SubSurf(1).ChordMesh=[2];%nchord1, nchord2
AC(2).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%%Htail
AC(3).Config.Proot=[11.5 0 1.5];%m  [x  y  z]
AC(3).Config.Ptip=[11.5 3  1.5];%m   [x  y  z]
AC(3).Config.rChord=1;%m
AC(3).Config.tChord=1;%m
AC(3).Config.Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(3).Config.Airfoil='flat';% Airfoil file
AC(3).Config.Dihedral=0*pi/180;% dihedral angle
AC(3).Config.Dihedral_ref = 1;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(3).Config.Symmetry='no';% symmetrical wing

AC(3).SubSurf(1).Name='Htail';% wing surface
AC(3).SubSurf(1).SpanData=[0 1 2];% span percentage interval, nspan
AC(3).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(3).SubSurf(1).ChordMesh=[2];%nchord1, nchord2
AC(3).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%%
AcControl = [];
AcTrim = [];