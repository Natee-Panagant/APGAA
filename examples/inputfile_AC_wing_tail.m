function [AC,FC] = inputfile_AC_wing_tail

%% 

%%
% Flight conditions
FC.Qinf = [10 0 0];% Velocity for static loadings, m/s
FC.rho_air = 1.20;% kg/m^3
FC.M = 0.5;% Mach number
FC.rG = [1.00 0 0];% center of mass
FC.k = [0.001 0.6 1.4]; %Nastran reduce frequencies (omega*Uinf/semichord)


%%
% Aircraft configuration details

%%
% Wing
wingXpos = 3; %wing x position, meter
wingYpos = 0.6;%wing y position, meter
wingZpos = -0.4;%wing z position, meter

AC(1).Label = 'Wing';
AC(1).Proot = [wingXpos wingYpos wingZpos];%m
AC(1).Ptip = [wingXpos+0.4 wingYpos+10 wingZpos];%m
AC(1).rChord = 1.25;    %m
AC(1).tChord = 1.2;      %m
AC(1).nSpanPanel = 20;% no. of spanwise vortex panels
AC(1).nChordPanel = 5;% no. of chordwise vortex panels
AC(1).Dihedral = 5*pi/180;% dihedral angle
AC(1).dRefPt = [0 0 wingZpos];% % reference point for dihedral angle rotating
AC(1).Incidence = 3*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(1).iRefPt = [wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(1).Airfoil = 'flat';% Airfoil file
AC(1).Symmetry = 'yes';% symmetrical wing  

% Horizontal tail
HtailXpos = 8.5;%Horizontal tail x position, meter
HtailYpos = 0.5;%Horizontal tail y position, meter
HtailZpos = 0.25;%Horizontal tail z position, meter

AC(2).Label = 'Horizontal tail';
AC(2).Proot = [HtailXpos HtailYpos HtailZpos];%m
AC(2).Ptip = [HtailXpos+0.4 HtailYpos+3 HtailZpos];%m
AC(2).rChord = .75;%m
AC(2).tChord = 0.5;%m
AC(2).nSpanPanel = 4;% no. of spanwise vortex panels
AC(2).nChordPanel = 4;% no. of chordwise vortex panels
AC(2).Dihedral = -5*pi/180;% dihedral angle
AC(2).dRefPt = [0 0 HtailZpos];% % reference point for dihedral angle rotating
AC(2).Incidence = -2*pi/180;% incidence, radian
AC(2).iRefPt = [HtailXpos HtailYpos HtailZpos];% reference point for incidence angle rotating
AC(2).Airfoil = 'flat';% Airfoil file
AC(2).Symmetry = 'yes';% symmetrical wing

% Vertical tail
VtailXpos = 8.5;%Vertical tail x position, meter
VtailYpos = 0.5;%Vertical tail x position, meter
VtailZpos = 0.0;%Vertical tail z position, meter
AC(3).Label = 'Vertical tail';
AC(3).Proot = [VtailXpos VtailYpos VtailZpos];%m
AC(3).Ptip = [VtailXpos+0.4 VtailYpos+2.5 VtailZpos];%m
AC(3).rChord = .75;%m
AC(3).tChord = 0.5;%m
AC(3).nSpanPanel = 4;% no. of spanwise vortex panels
AC(3).nChordPanel = 4;% no. of chordwise vortex panels
AC(3).Dihedral = 90*pi/180;% dihedral angle
AC(3).dRefPt = [0 0 VtailZpos];% % reference point for dihedral angle rotating
AC(3).Incidence = 0*pi/180;% incidence, radian
AC(3).iRefPt = [VtailXpos VtailYpos VtailZpos];% reference point for incidence angle rotating
AC(3).Airfoil = 'flat';% Airfoil file
AC(3).Symmetry = 'no';% Non-symmetric wing



end
