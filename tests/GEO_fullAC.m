function [AC, AcControl, AcTrim]=GEO_fullAC(state)
%% Wing section 1 horizontal
AC(1).Label='Wing';
AC(1).State.Qinf = state.Qinf;% Velocity for static loadings, m/s
AC(1).State.rho_air=state.rho_air;%kg/m^3
AC(1).State.M=state.M;% Mach number
AC(1).State.dt=0.5;% time interval
AC(1).State.rG=state.CG;% center of mass
AC(1).State.k=state.k; %Nastran reduce frequencies (omega*Uinf/semichord)

AC(1).Config.Proot=[21 2 0];%m  [x  y  z]
AC(1).Config.Ptip=[17 15 0];%m   [x  y  z]
AC(1).Config.rChord=5;%m
AC(1).Config.tChord=2;%m
AC(1).Config.Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).Config.Airfoil='flat';% Airfoil file
AC(1).Config.Dihedral=0*pi/180;% dihedral angle
AC(1).Config.Dihedral_ref = 1;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(1).Config.Symmetry='yes';% symmetrical wing

% Aero Analysis (Gen Panel)                                                       
AC(1).SubSurf(1).Name='Wing';% wing surface
AC(1).SubSurf(1).SpanData=[0 1 10];% span percentage interval, nspan
AC(1).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(1).SubSurf(1).ChordMesh=[4];%nchord1, nchord2
AC(1).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%% Vertical tail
AC(2).Label='Vertical tail';
AC(2).State = AC(1).State;
AC(2).Config.Proot=[32 0 2];%m (*Vertical tail must input in x-y plane then rotate with dihedral angle)
AC(2).Config.Ptip=[37 0 7];%m  (*Vertical tail must input in x-y plane then rotate with dihedral angle)
AC(2).Config.rChord=6;%m
AC(2).Config.tChord=3;%m
AC(2).Config.Incidence=0*pi/180;% incidence, radian
AC(2).Config.Airfoil='flat';% Airfoil file
AC(2).Config.Dihedral=0*pi/180;% dihedral angle
AC(2).Config.Dihedral_ref = 2;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(2).Config.Symmetry='no';% symmetrical wing

% Aero Analysis (Gen Panel)   
AC(2).SubSurf(1).Name='Vertical Tail';% aileron surface
AC(2).SubSurf(1).SpanData=[0 1 8];% span percentage interval, nspan
AC(2).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(2).SubSurf(1).ChordMesh=[4];%nchord1, nchord2
AC(2).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%% Horizontal tail
AC(3).Label='Horizontal tail';
AC(3).State = AC(1).State;
AC(3).Config.Proot=[37 0 7];%m
AC(3).Config.Ptip=[40 7 7];%m
AC(3).Config.rChord=3;%m
AC(3).Config.tChord=2;%m
AC(3).Config.Incidence=0*pi/180;% incidence, radian
AC(3).Config.Airfoil='flat';% Airfoil file
AC(3).Config.Dihedral=0*pi/180;% dihedral angle
AC(3).Config.Dihedral_ref = 1;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(3).Config.Symmetry='yes';% symmetrical wing

% Aero Analysis (Gen Panel)   
AC(3).SubSurf(1).Name='Horizontal tail';% aileron surface
AC(3).SubSurf(1).SpanData=[0 1 8];% span percentage interval, nspan
AC(3).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(3).SubSurf(1).ChordMesh=[4];%nchord1, nchord2
AC(3).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%% Fuselage
AC(4).Label='Fuselage';
AC(4).State = AC(1).State;
AC(4).Config.Proot=[0 0 0];%m
AC(4).Config.Ptip=[0 2 0];%m
AC(4).Config.rChord=40;%m
AC(4).Config.tChord=40;%m
AC(4).Config.Incidence=0*pi/180;% incidence, radian
AC(4).Config.Airfoil='flat';% Airfoil file
AC(4).Config.Dihedral=0*pi/180;% dihedral angle
AC(4).Config.Dihedral_ref = 1;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(4).Config.Symmetry='yes';% symmetrical wing

% Aero Analysis (Gen Panel)   
AC(4).SubSurf(1).Name='Fuselage';% aileron surface
AC(4).SubSurf(1).SpanData=[0 1 1];% span percentage interval, nspan
AC(4).SubSurf(1).ChordIntv=[0 1];%chord percentage interval
AC(4).SubSurf(1).ChordMesh=[20];%nchord1, nchord2
AC(4).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

AcControl = [];
AcTrim = [];