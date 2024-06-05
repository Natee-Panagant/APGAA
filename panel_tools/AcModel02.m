function AC=AcModel02
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
wingZpos=-0.2;%wing z position, meter
AC(2).Label='Wing';
AC(2).Proot=[wingXpos wingYpos wingZpos];%m
AC(2).Ptip=[wingXpos+0.4 wingYpos+5 wingZpos];%m
AC(2).rChord=2.25;%m
AC(2).tChord=1.2;%m
AC(2).nSpanPanel=12;% no. of spanwise vortex panels
AC(2).nChordPanel=5;% no. of chordwise vortex panels
AC(2).Incidence=5*pi/180;% incidence, radian and rotating point (Proot of surface 1)
AC(2).iRefPt=[wingXpos wingYpos wingZpos];% reference point for incidence angle rotating
AC(2).Airfoil='N64210.dat';% Airfoil file
AC(2).Dihedral=5*pi/180;% dihedral angle
AC(2).dRefPt=[0 0 wingZpos];% % reference point for dihedral angle rotating 
AC(2).Symmetry='yes';% symmetrical wing  
AC(2).nSpan=10;% no. of spanwise sections
AC(2).nChord=15;% no. of chordwise sections

AC(2).SubSurf(1).Name='Flap';% flap surface
AC(2).SubSurf(1).SpanData=[0.2 0.5];% span percentage interval, nspan
% ChordIntv=[main wing leading edge, sub-surface leading edge, wing trailing edge];
AC(2).SubSurf(1).ChordIntv=[0.0 0.7 1.0];%chord percentage interval
AC(2).SubSurf(1).ChordMesh=5;%nchord1, nchord2

AC(2).SubSurf(2).Name='Aileron';% aileron surface
AC(2).SubSurf(2).SpanData=[0.55 1.0];% span percentage interval, nspan
AC(2).SubSurf(2).ChordIntv=[0.0 0.7 1.0];%chord percentage interval
AC(2).SubSurf(2).ChordMesh=4;%nchord1, nchord2

% Horizontal tail
HtailXpos=8.5;%Horizontal tail x position, meter
HtailYpos=0.5;%Horizontal tail y position, meter
HtailZpos=0.25;%Horizontal tail z position, meter
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
AC(3).Dihedral=35*pi/180;% dihedral angle
AC(3).dRefPt=[0 0 HtailZpos];% % reference point for dihedral angle rotating
AC(3).Symmetry='yes';% symmetrical wing
AC(3).nSpan=10;% no. of spanwise sections
AC(3).nChord=15;% no. of chordwise sections

AC(3).SubSurf(1).Name='Elevator';% aileron surface
AC(3).SubSurf(1).SpanData=[0.0 1.0 4];% span percentage interval, nspan
AC(3).SubSurf(1).ChordIntv=[0.0 0.6 1.0];%chord percentage interval
AC(3).SubSurf(1).ChordMesh=4;%nchord for sub-surface


