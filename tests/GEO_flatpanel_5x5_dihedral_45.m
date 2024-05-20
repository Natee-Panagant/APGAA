function [AC, AcControl, AcTrim]=GEO_flatpanel_5x5_dihedral_45(state)
%% Wing section 1 horizontal
AC(1).Label='Wing';
AC(1).State.Qinf = state.Qinf;% Velocity for static loadings, m/s
AC(1).State.rho_air=state.rho_air;%kg/m^3
AC(1).State.M=state.M;% Mach number
AC(1).State.dt=0.5;% time interval
AC(1).State.rG=state.CG;% center of mass
AC(1).State.k=state.k; %Nastran reduce frequencies (omega*Uinf/semichord)

AC(1).Config.Proot=[0 0 0];%m  [x  y  z]
AC(1).Config.Ptip=[0 4 0];%m   [x  y  z]
AC(1).Config.rChord=1;%m
AC(1).Config.tChord=1;%m
AC(1).Config.Incidence=0*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).Config.Airfoil='flat';% Airfoil file
AC(1).Config.Dihedral=45*pi/180;% dihedral angle
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