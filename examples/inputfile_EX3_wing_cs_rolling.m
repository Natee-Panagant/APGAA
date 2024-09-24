function [AC,FC] = inputfile_EX3_wing_cs_rolling
% Flight conditions
FC.Qinf = [10 , 0 , 0];% Velocity for static loadings, m/s
FC.rho_air = 1.20;% kg/m^3
FC.M = 0.5;% Mach number
FC.rG = [3.25 , 0 , 0];% center of mass
FC.k = [0.001 , 0.6 , 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)

%% AIRCRAFT CONFIGURATION DETAILS

%%%%%%%%
% Wing %
%%%%%%%%
AC(1).Label = 'Wing';
AC(1).Proot = [3 , 0 , -0.4];%(m) [x,y,z] Position of point at root 
AC(1).Ptip = [3.4 , 6 , -0.4];%% (m) [x,y,z] Position of point at tip 
AC(1).rChord = 1.2;% (m) Root chord
AC(1).tChord = 0.8;% (m) Tip chord
AC(1).nSpanPanel = 20;% (-) Number of vortex panels in spanwise
AC(1).nChordPanel = 5;% (-) Number of vortex panels in chordwise

AC(1).Dihedral = 2*pi/180;% (rad) Dihedral angle
AC(1).dRefPt = [0 , 0 , -0.4];% (m) [x,y,z] Reference point for dihedral angle for rotating
AC(1).Incidence = 2*pi/180;% (rad) Incidence angle
AC(1).iRefPt = [3 , 0 , -0.4];% (m) [x,y,z] Reference point for incidence angle for rotating
AC(1).Airfoil = 'flat';% (-) [string] Airfoil profile name
AC(1).Symmetry = 'yes';% (yes/no) Symmetrical wing behavior

% Flap
AC(1).SubSurf(1).Name = 'Flap';% Label of sub-surf: Flap
AC(1).SubSurf(1).SpanData = [0.1 , 0.5 , 1];% (%) [rspan1 rspan2 ... 1] Percentage interval spanwise on surf 1
AC(1).SubSurf(1).ChordIntv = [0.0 , 0.7 , 1.0];% (%) [rchord1 rchord2 ... 1] Percentage interval chordwise on surf 1
AC(1).SubSurf(1).ChordMesh = 1;% (-) Number of chordwise of sub-surf no. 1
AC(1).SubSurf(1).RotAngle = 0*pi/180;% (rad) Rotation angle of control surface at sub-surf no. 1
% Aileron
AC(1).SubSurf(2).Name = 'Aileron';% Label of sub-surf: Aileron
AC(1).SubSurf(2).SpanData = [0.5 , 1.0 , 1];% (%) [rspan1 rspan2 ... 1] Percentage interval spanwise on surf 1
AC(1).SubSurf(2).ChordIntv = [0.0 , 0.7 , 1.0];% (%) [rchord1 rchord2 ... 1] Percentage interval chordwise on surf 1
AC(1).SubSurf(2).ChordMesh = 1;% (-) Number of chordwise of sub-surf no. 2
AC(1).SubSurf(2).RotAngle = 5*pi/180;% (rad) Rotation angle of control surface at sub-surf no. 2

end