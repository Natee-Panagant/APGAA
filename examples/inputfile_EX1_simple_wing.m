function [AC,FC] = inputfile_EX1_simple_wing
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

end