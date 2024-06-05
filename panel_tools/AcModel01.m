function [AC,FC]=AcModel01
% Flight conditions
FC.Qinf=[10 0 0];% Velocity for static loadings, m/s
FC.rho_air=1.20;%kg/m^3
FC.Mach=0.5;% Mach number
FC.dt=0.1;% time interval
FC.nWakePanel=4;% no. of chordwise wake panels
FC.rG=[1.00 0 0];% center of mass

% Aircraft configuration details
% Fuselage
AC(1).Label='Fuselage';
AC(1).Pstart=[0 0 0];% Point for fuselage nose
AC(1).Length=10;% Fuselage length
AC(1).xStations=linspace(0,1,10);%Normalised axial positions of predefined cross-sections 
AC(1).a=0.6*[0.01 0.75 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.5];% y-axis widths for super-ellipses
AC(1).bU=0.5*[0.1 1.15 2.00 1.45 1.25 1.1 1.0 1.0 1.0 0.5];% upper surface heights for super-ellipses
AC(1).bL=-0.5*[0.1 0.75 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.5];% lower surface heights for super-ellipses
AC(1).kU=4.55*ones(1,10);% upper surface powers for super-ellipses
AC(1).kL=4.55*ones(1,10);% lower surface powers for super-ellipses
AC(1).RadialResol=16;% Theta-axis resolution
AC(1).AxialResol=20;% Axial resolution

% Wing
wingXpos=3;%wing x position, meter
wingYpos=0.6;%wing y position, meter
wingZpos=-0.4;%wing z position, meter
AC(2).Label='Wing';
AC(2).Proot=[wingXpos wingYpos wingZpos];%m
AC(2).Ptip=[wingXpos+0.4 wingYpos+10 wingZpos];%m
AC(2).rChord=1.25;%m
AC(2).tChord=1.2;%m
AC(2).nSpanPanel=5;% no. of spanwise vortex panels
AC(2).nChordPanel=4;% no. of chordwise vortex panels
AC(2).Incidence=3*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(2).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(2).Airfoil='N64210.dat';% Airfoil file
AC(2).Dihedral=5*pi/180;% dihedral angle
AC(2).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(2).Symmetry='yes';% symmetrical wing  
AC(2).nSpan=10;% no. of spanwise sections for wing surfaces
AC(2).nChord=15;% no. of chordwise sections for wing surfaces

AC(2).SubSurf(1).Name='Flap';% flap surface
AC(2).SubSurf(1).SpanData=[0.1 0.6 4];% span percentage interval, nspan
AC(2).SubSurf(1).ChordIntv=[0.0 0.7 1.0];%chord percentage interval
AC(2).SubSurf(1).ChordMesh=4;%nchord1, nchord2
AC(2).SubSurf(2).Name='Aileron';% aileron surface
AC(2).SubSurf(2).SpanData=[0.6 1.0];% span percentage interval, nspan
AC(2).SubSurf(2).ChordIntv=[0.0 0.7 1.0];%chord percentage interval
AC(2).SubSurf(2).ChordMesh=5;%nchord1, nchord2

% Horizontal tail
HtailXpos=8.5;%Horizontal tail x position, meter
HtailYpos=0.5;%Horizontal tail y position, meter
HtailZpos=0.0;%Horizontal tail z position, meter
AC(3).Label='Horizontal tail';
AC(3).Proot=[HtailXpos HtailYpos HtailZpos];%m
AC(3).Ptip=[HtailXpos+0.4 HtailYpos+3 HtailZpos];%m
AC(3).rChord=.75;%m
AC(3).tChord=0.5;%m
AC(3).nSpanPanel=5;% no. of spanwise vortex panels
AC(3).nChordPanel=4;% no. of chordwise vortex panels
AC(3).Incidence=-2*pi/180;% incidence, radian
AC(3).iRefPt=[HtailXpos HtailYpos HtailZpos];% reference point for incidence angle rotating
AC(3).Airfoil='naca0012.DAT';% Airfoil file
AC(3).Dihedral=0*pi/180;% dihedral angle
AC(3).dRefPt=[0 0 HtailZpos];% % reference point for dihedral angle rotating
AC(3).Symmetry='yes';% symmetrical wing
AC(3).nSpan=10;% no. of spanwise sections for wing surfaces
AC(3).nChord=15;% no. of chordwise sections for wing surfaces

AC(3).SubSurf(1).Name='Elevator';% aileron surface
AC(3).SubSurf(1).SpanData=[0.0 1.0];% span percentage interval, nspan
AC(3).SubSurf(1).ChordIntv=[0.0 0.6 1.0];%chord percentage interval
AC(3).SubSurf(1).ChordMesh=3;%nchord1, nchord2

% Vertical tail
VtailXpos=8.5;%Vertical tail x position, meter
VtailYpos=0.5;%Vertical tail x position, meter
VtailZpos=0.0;%Vertical tail z position, meter
AC(4).Label='Vertical tail';
AC(4).Proot=[VtailXpos VtailYpos VtailZpos];%m
AC(4).Ptip=[VtailXpos+0.4 VtailYpos+2.5 VtailZpos];%m
AC(4).rChord=.75;%m
AC(4).tChord=0.5;%m
AC(4).nSpanPanel=4;% no. of spanwise vortex panels
AC(4).nChordPanel=4;% no. of chordwise vortex panels
AC(4).Incidence=0*pi/180;% incidence, radian
AC(4).iRefPt=[VtailXpos VtailYpos VtailZpos];% reference point for incidence angle rotating
AC(4).Airfoil='naca0012.DAT';% Airfoil file
AC(4).Dihedral=90*pi/180;% dihedral angle
AC(4).dRefPt=[0 0 VtailZpos];% % reference point for dihedral angle rotating
AC(4).Symmetry='no';% Non-symmetric wing
AC(4).nSpan=10;% no. of spanwise sections for wing surfaces
AC(4).nChord=15;% no. of chordwise sections for wing surfaces

AC(4).SubSurf(1).Name='Rudder';% aileron surface
AC(4).SubSurf(1).SpanData=[0.0 1.0];% span percentage interval, nspan
AC(4).SubSurf(1).ChordIntv=[0.0 0.6 1.0];%chord percentage interval
AC(4).SubSurf(1).ChordMesh=3;%



