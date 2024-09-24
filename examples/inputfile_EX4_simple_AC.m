function [AC,FC] = inputfile_EX4_simple_AC
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

%%%%%%%%%%%%%%%%%%%
% Horizontal tail %
%%%%%%%%%%%%%%%%%%%
AC(2).Label = 'Horizontal tail';
AC(2).Proot = [8.5 , 0 , 0.25];%m
AC(2).Ptip = [9 , 2.5 , 0.25];%m
AC(2).rChord = 1;%m
AC(2).tChord = 0.5;%m
AC(2).nSpanPanel = 4;% no. of spanwise vortex panels
AC(2).nChordPanel = 4;% no. of chordwise vortex panels
AC(2).Dihedral = 0*pi/180;% dihedral angle
AC(2).dRefPt = [0 0 0.25];% % reference point for dihedral angle rotating
AC(2).Incidence = -1*pi/180;% incidence, radian
AC(2).iRefPt = [8.5 , 0 , 0.25];% reference point for incidence angle rotating
AC(2).Airfoil = 'flat';% Airfoil file
AC(2).Symmetry = 'yes';% symmetrical wing

%%%%%%%%%%%%%%%%%
% Vertical tail %
%%%%%%%%%%%%%%%%%
AC(3).Label = 'Vertical tail';
AC(3).Proot = [8.5 , 0 , 0.25];%m
AC(3).Ptip = [9 , 0 , 2.75];%m
AC(3).rChord = 1;%m
AC(3).tChord = 0.5;%m
AC(3).nSpanPanel = 4;% no. of spanwise vortex panels
AC(3).nChordPanel = 4;% no. of chordwise vortex panels
AC(3).Dihedral = 0*pi/180;% dihedral angle
AC(3).dRefPt = [0 , 0 , 0.25];% % reference point for dihedral angle rotating
AC(3).Incidence = 0*pi/180;% incidence, radian
AC(3).iRefPt = [8.5 , 0 , 0.25];% reference point for incidence angle rotating
AC(3).Airfoil = 'flat';% Airfoil file
AC(3).Symmetry = 'no';% Non-symmetric wing

end