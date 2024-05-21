function AC = ex_simple_wing
%% Wing
AC(1).Label='Wing';
AC(1).State.Qinf = [100*cosd(5) 100*sind(5) 0];% Velocity [Vx Vy Vz] for static loadings, m/s
AC(1).State.rho_air = 1.225;%kg/m^3
AC(1).State.M = 0.8;% Mach number
AC(1).State.k = [0.001 0.6 1.4];; %Nastran reduce frequencies (omega*Uinf/semichord)

AC(1).Config.Proot = [4.5 0 0];%m  [x  y  z]
AC(1).Config.Ptip = [4.5 3 0];%m   [x  y  z]
AC(1).Config.rChord = 1.5;%m
AC(1).Config.tChord = 1.5;%m
AC(1).Config.Incidence = 0*pi/180;% incidence (radian)
AC(1).Config.Airfoil = 'flat';% Airfoil file
AC(1).Config.Dihedral = 0*pi/180;% dihedral angle (radian)
AC(1).Config.Dihedral_ref = 1;% = 1 -> rotate around x-axis, = 2 -> rotate around Proot
AC(1).Config.Symmetry = 'no';% symmetrical wing

%Lifting Surface
AC(1).SubSurf(1).Name = 'Wing';% wing surface
AC(1).SubSurf(1).SpanData = [0 1 2];% span percentage interval, nspan
AC(1).SubSurf(1).ChordIntv = [0 0.8];%chord percentage interval
AC(1).SubSurf(1).ChordMesh = [2];%nchord1, nchord2
AC(1).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%Control Surface
AC(1).SubSurf(2).Name = 'Flap';% wing surface
AC(1).SubSurf(2).SpanData = [0 1 2];% span percentage interval, nspan
AC(1).SubSurf(2).ChordIntv = [0.8 1];%chord percentage interval
AC(1).SubSurf(2).ChordMesh = [2];%nchord1, nchord2
AC(1).SubSurf(2).type = 2; % 1 = non-control surface, 2 = control surface
AC(1).SubSurf(2).adjust_angle = 5*pi/180;%adjustable angle for control surface
