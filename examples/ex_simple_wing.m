function [AC,State] = ex_simple_wing
State.Qinf = [274.4 0 0];% Velocity [Vx Vy Vz] for static loadings, m/s
State.rho_air = 1.225;%kg/m^3
State.M = 0.8;% Mach number
State.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)

% Wing
AC(1).Label='Wing';
AC(1).Config.Proot = [4.5 0 0];%m  [x  y  z]
AC(1).Config.Ptip = [5.5 3 0];%m   [x  y  z]
AC(1).Config.rChord = 1.5;%m
AC(1).Config.tChord = 1;%m
AC(1).Config.Incidence = 5*pi/180;% incidence (radian)
AC(1).Config.Airfoil = 'flat';% Airfoil file
AC(1).Config.Dihedral = 15*pi/180;% dihedral angle (radian)
AC(1).Config.Dihedral_ref = 1;% = 1 -> rotate around origin, = 2 -> rotate around Proot
AC(1).Config.Symmetry = 'yes';% symmetrical wing

%Lifting Surface
AC(1).SubSurf(1).Name = 'Wing';% wing surface
AC(1).SubSurf(1).SpanData = [0 1 5];% span percentage interval, nspan
AC(1).SubSurf(1).ChordIntv = [0 0.8];%chord percentage interval
AC(1).SubSurf(1).ChordMesh = [4];%nchord1, nchord2
AC(1).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%Control Surface
AC(1).SubSurf(2).Name = 'Aileron';% wing surface
AC(1).SubSurf(2).SpanData = [0 1 5];% span percentage interval, nspan
AC(1).SubSurf(2).ChordIntv = [0.8 1];%chord percentage interval
AC(1).SubSurf(2).ChordMesh = [1];%nchord1, nchord2
AC(1).SubSurf(2).type = 2; % 1 = non-control surface, 2 = control surface
AC(1).SubSurf(2).adjust_angle = 10*pi/180;%adjustable angle for control surface

% Vtail
AC(2).Label='Vtail';
AC(2).Config.Proot = [7 0 0];%m  [x  y  z]
AC(2).Config.Ptip = [8 0 1];%m   [x  y  z]
AC(2).Config.rChord = 1;%m
AC(2).Config.tChord = 0.5;%m
AC(2).Config.Incidence = 0*pi/180;% incidence (radian)
AC(2).Config.Airfoil = 'flat';% Airfoil file
AC(2).Config.Dihedral = 0*pi/180;% dihedral angle (radian)
AC(2).Config.Dihedral_ref = 1;% = 1 -> rotate around origin, = 2 -> rotate around Proot
AC(2).Config.Symmetry = 'no';% symmetrical wing

%Vertical stabilizer
AC(2).SubSurf(1).Name = 'Vtail';% wing surface
AC(2).SubSurf(1).SpanData = [0 1 5];% span percentage interval, nspan
AC(2).SubSurf(1).ChordIntv = [0 0.8];%chord percentage interval
AC(2).SubSurf(1).ChordMesh = [4];%nchord1, nchord2
AC(2).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%Control Surface
AC(2).SubSurf(2).Name = 'Rudder';% wing surface
AC(2).SubSurf(2).SpanData = [0 1 5];% span percentage interval, nspan
AC(2).SubSurf(2).ChordIntv = [0.8 1];%chord percentage interval
AC(2).SubSurf(2).ChordMesh = [1];%nchord1, nchord2
AC(2).SubSurf(2).type = 2; % 1 = non-control surface, 2 = control surface
AC(2).SubSurf(2).adjust_angle = 15*pi/180;%adjustable angle for control surface

% Htail
AC(3).Label='Htail';
AC(3).Config.Proot = [7 0 0];%m  [x  y  z]
AC(3).Config.Ptip = [8 1 0];%m   [x  y  z]
AC(3).Config.rChord = 0.8;%m
AC(3).Config.tChord = 0.5;%m
AC(3).Config.Incidence = 0*pi/180;% incidence (radian)
AC(3).Config.Airfoil = 'flat';% Airfoil file
AC(3).Config.Dihedral = 10*pi/180;% dihedral angle (radian)
AC(3).Config.Dihedral_ref = 2;% = 1 -> rotate around origin, = 2 -> rotate around Proot
AC(3).Config.Symmetry = 'yes';% symmetrical wing

%Vertical stabilizer
AC(3).SubSurf(1).Name = 'Htail';% wing surface
AC(3).SubSurf(1).SpanData = [0 1 5];% span percentage interval, nspan
AC(3).SubSurf(1).ChordIntv = [0 0.8];%chord percentage interval
AC(3).SubSurf(1).ChordMesh = [4];%nchord1, nchord2
AC(3).SubSurf(1).type = 1; % 1 = non-control surface, 2 = control surface

%Control Surface
AC(3).SubSurf(2).Name = 'Elevator';% wing surface
AC(3).SubSurf(2).SpanData = [0 1 5];% span percentage interval, nspan
AC(3).SubSurf(2).ChordIntv = [0.8 1];%chord percentage interval
AC(3).SubSurf(2).ChordMesh = [1];%nchord1, nchord2
AC(3).SubSurf(2).type = 2; % 1 = non-control surface, 2 = control surface
AC(3).SubSurf(2).adjust_angle = 15*pi/180;%adjustable angle for control surface