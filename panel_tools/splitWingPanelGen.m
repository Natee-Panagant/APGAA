function WgOut = splitWingPanelGen(Wg,iSurf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating panels for a surface with sub-surfaces (control surfaces) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transform vertical Ptip input to horizontal
if Wg.Proot(2) == Wg.Ptip(2)
    Wg.Ptip(2) = Wg.Proot(2)+(Wg.Ptip(3)-Wg.Proot(3));
    Wg.Ptip(3) = Wg.Proot(3);
    Wg.Dihedral = 90*pi/180;
    Wg.dRefPt = Wg.Proot;
end

% Load airfoil data
if strcmp(Wg.Airfoil,'flat')
    Airfoil = [linspace(0,1,20)',zeros(20,1)];
else
    Airfoil = load(Wg.Airfoil);
end
% Split Upper and Lower surfaces
if Airfoil(1,1)>0
    UpperX = Airfoil(2:Airfoil(1,1)+1,1);
    UpperZ = Airfoil(2:Airfoil(1,1)+1,2);
    LowerZ = Airfoil(Airfoil(1,1)+2:end,2);
else
    nAfoil = size(Airfoil,1);
    UpperX = Airfoil(1:nAfoil/2,1);
    UpperZ = Airfoil(1:nAfoil/2,2);
    LowerZ = Airfoil(nAfoil/2+1:end,2);
end

% Split the main surface for generating control surfaces
chordRatio = Wg.SubSurf(iSurf).ChordIntv(2);
nChordM = max(1,round(chordRatio*Wg.nChordPanel));%M for main wing
spStart = Wg.SubSurf(iSurf).SpanData(1);
spEnd = Wg.SubSurf(iSurf).SpanData(2);

[UpperX,nsort] = sort(UpperX);
UpperZ = UpperZ(nsort);
LowerZ = LowerZ(nsort);
CamberZ = 0.5*(UpperZ+LowerZ);

% Convert airfoil data from [0,1] to real scale
maxX = max(UpperX); minX = min(UpperX);
dX = maxX-minX;
delXm = 0.25*(maxX-minX)/Wg.nChordPanel;
delXs = 0.25*(1-chordRatio)*(maxX-minX)/Wg.SubSurf(iSurf).ChordMesh;
xmWing = linspace(minX,minX+chordRatio*dX,nChordM+1);
xsWing = linspace(minX+chordRatio*dX,maxX,Wg.SubSurf(iSurf).ChordMesh+1);
xWing = [xmWing xsWing(2:end)];
zWing = interp1(UpperX,CamberZ,xWing,'linear','extrap');
xmPanel = linspace(minX+delXm,(minX+chordRatio*dX)+delXm,nChordM+1);
xsPanel = linspace((minX+chordRatio*dX)+delXs,maxX+delXs,Wg.SubSurf(iSurf).ChordMesh+1);
xPanel = [xmPanel(1:end-1) xsPanel(1:end-1) maxX+delXm];
zPanel = interp1(UpperX,CamberZ,xPanel,'linear','extrap');

Wgm.nChordM = nChordM;
Wgm.xWing = (xWing-min(xWing))/(max(xWing)-min(xWing));
Wgm.zWing = zWing/(max(xWing)-min(xWing));
Wgm.xPanel = (xPanel-min(xPanel))/(max(xPanel)-min(xPanel));
Wgm.zPanel = zPanel/(max(xPanel)-min(xPanel));
Wgm.Proot = [Wg.Proot(1)+spStart*(Wg.Ptip(1)-Wg.Proot(1)),...
             Wg.Proot(2)+spStart*(Wg.Ptip(2)-Wg.Proot(2)),...
             Wg.Proot(3)+spStart*(Wg.Ptip(3)-Wg.Proot(3))];
Wgm.Ptip = [Wg.Proot(1)+spEnd*(Wg.Ptip(1)-Wg.Proot(1)),...
            Wg.Proot(2)+spEnd*(Wg.Ptip(2)-Wg.Proot(2)),...
            Wg.Proot(3)+spEnd*(Wg.Ptip(3)-Wg.Proot(3))];
Wgm.rChord = (Wg.rChord+spStart*(Wg.tChord-Wg.rChord));
Wgm.tChord = (Wg.rChord+spEnd*(Wg.tChord-Wg.rChord));
Wgm.chordRatio = chordRatio;
Wgm.dXroot = 0.25*Wgm.rChord/Wg.nChordPanel;% leading root vortex ring position
Wgm.dXtip = 0.25*Wgm.tChord/Wg.nChordPanel;
Wgm.nChordPanel = Wg.SubSurf(iSurf).ChordMesh+nChordM;
Wgm.nSpanPanel = max(1,round(abs(spEnd-spStart)*Wg.nSpanPanel));
Wgm.Incidence = Wg.Incidence;
Wgm.iRefPt = Wg.iRefPt;
Wgm.Dihedral = Wg.Dihedral;
Wgm.dRefPt = Wg.dRefPt;
% Generate surface grids
[Xw,Yw,Zw,Xp,Yp,Zp,pHinge,rotAxis] = GenWingGrid(Wgm);%w for wings, p for panels
% Generate Hinge and rotation axis
pnNum = reshape(1:Wgm.nChordPanel*Wgm.nSpanPanel,Wgm.nChordPanel,Wgm.nSpanPanel);
subSurfNum = reshape(pnNum(nChordM+1:Wgm.nChordPanel,:),1,...
    Wg.SubSurf(iSurf).ChordMesh*Wgm.nSpanPanel);
WgOut.RpanelNum = subSurfNum;%sub-surface (control surface) panel numbers
WgOut.HingePt_R = pHinge;% sub-surface hinge position
WgOut.RotAxis_R = rotAxis;% sub-surface rotating axis

WgOut.Right.Xw = Xw; WgOut.Right.Yw = Yw; WgOut.Right.Zw = Zw;
WgOut.Right.Xp = Xp; WgOut.Right.Yp = Yp; WgOut.Right.Zp = Zp;
if strcmp(Wg.Symmetry,'yes')
    WgOut.Left.Xw = Xw; WgOut.Left.Yw = -Yw; WgOut.Left.Zw = Zw;
    WgOut.Left.Xp = Xp; WgOut.Left.Yp = -Yp; WgOut.Left.Zp = Zp;
    WgOut.LpanelNum = Wgm.nChordPanel*Wgm.nSpanPanel+subSurfNum;
    WgOut.HingePt_L = [pHinge(1) -pHinge(2) pHinge(3)];% sub-surface hinge position
    WgOut.RotAxis_L = [rotAxis(1) -rotAxis(2) rotAxis(3)];% sub-surface rotating axis
end

switch Wg.SubSurf(iSurf).Name
    case 'Flap'
        if WgOut.RotAxis_R(2)>0 %Rotation Axis of the Right-side Flap pointing inward (To the Left)
            WgOut.HingePt_R = WgOut.HingePt_R+WgOut.RotAxis_R;
            WgOut.RotAxis_R = -WgOut.RotAxis_R;
        end
        if strcmp(Wg.Symmetry,'yes')
            if WgOut.RotAxis_L(2)>0 %Rotation Axis of the Right-side Flap pointing inward (To the Left)
                WgOut.HingePt_L = WgOut.HingePt_L+WgOut.RotAxis_L;
                WgOut.RotAxis_L = -WgOut.RotAxis_L;%Rotation Axis of the Left-side Flap pointing outward (flip + inverse Y)
            end
        end
    case 'Aileron'
        if WgOut.RotAxis_R(2)<0 %Rotation Axis of the Right-side wing pointing outward (To the Right)
            WgOut.HingePt_R = WgOut.HingePt_R+WgOut.RotAxis_R;
            WgOut.RotAxis_R = -WgOut.RotAxis_R;
        end
        if strcmp(Wg.Symmetry,'yes')
            if WgOut.RotAxis_L(2)>0 %Rotation Axis of the Left-side wing pointing outward (inverse Y) (To the Left)
                WgOut.HingePt_L = WgOut.HingePt_L+WgOut.RotAxis_L;
                WgOut.RotAxis_L = -WgOut.RotAxis_L;
            end
        end
    case 'Elevator'
        if WgOut.RotAxis_R(2)>0 %Rotation Axis of the Right-side elevator pointing inward (To the Left)
            WgOut.HingePt_R = WgOut.HingePt_R+WgOut.RotAxis_R;
            WgOut.RotAxis_R = -WgOut.RotAxis_R;
        end
        if strcmp(Wg.Symmetry,'yes')
            if WgOut.RotAxis_L(2)>0 %Rotation Axis of the Left-side elevator pointing outward (flip + inverse Y) (To the Left)
                WgOut.HingePt_L = WgOut.HingePt_L+WgOut.RotAxis_L;
                WgOut.RotAxis_L = -WgOut.RotAxis_L;
            end
        end
    case 'Rudder'
        if WgOut.RotAxis_R(3)<0 %Rotation Axis of the rudder is pointing upward
            WgOut.HingePt_R = WgOut.HingePt_R+WgOut.RotAxis_R;
            WgOut.RotAxis_R = -WgOut.RotAxis_R;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Sub-Functions %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xw,Yw,Zw,Xp,Yp,Zp,pHinge,rotAxis] = GenWingGrid(Wg)
ProotW = Wg.Proot;
ProotP = [Wg.Proot(1)+Wg.dXroot Wg.Proot(2:3)];
PtipW = Wg.Ptip;
PtipP = [Wg.Ptip(1)+Wg.dXtip Wg.Ptip(2:3)];

Wg.nChordPanel = Wg.nChordPanel+1;% nChord intervals with nChord+1 nodes
Wg.nSpanPanel = Wg.nSpanPanel+1;% nSpan intervals with nSpan+1 nodes
x0W = linspace(ProotW(1),PtipW(1),Wg.nSpanPanel);
x0P = linspace(ProotP(1),PtipP(1),Wg.nSpanPanel);
yi = linspace(Wg.Proot(2),Wg.Ptip(2),Wg.nSpanPanel);
chDist = linspace(Wg.rChord,Wg.tChord,Wg.nSpanPanel);
for i = 1:Wg.nChordPanel
    Yw(i,:) = yi;
    Xw(i,:) = x0W+Wg.xWing(i)*chDist;
    Zw(i,:) = ProotW(3)+Wg.zWing(i)*chDist;
    Yp(i,:) = Yw(i,:);
    Xp(i,:) = x0P+Wg.xPanel(i)*chDist;
    Zp(i,:) = ProotP(3)+Wg.zPanel(i)*chDist;
end

% Move Hinge from structural grid to aerodynamic grid
pHingeR = [Xp(Wg.nChordM+1,1);Yp(Wg.nChordM+1,1);Zp(Wg.nChordM+1,1)];
pHingeT = [Xp(Wg.nChordM+1,end);Yp(Wg.nChordM+1,end);Zp(Wg.nChordM+1,end)];
% Rotate Panels with dihedral angle
if Wg.Dihedral ~= 0
    [Xw,Yw,Zw] = AxisRotating(Xw,Yw,Zw,Wg.Dihedral,[1;0;0],Wg.dRefPt');
    [Xp,Yp,Zp] = AxisRotating(Xp,Yp,Zp,Wg.Dihedral,[1;0;0],Wg.dRefPt');
    [pHingeR(1),pHingeR(2),pHingeR(3)] = AxisRotating(pHingeR(1),pHingeR(2),...
        pHingeR(3),Wg.Dihedral,[1;0;0],Wg.dRefPt');
    [pHingeT(1),pHingeT(2),pHingeT(3)] = AxisRotating(pHingeT(1),pHingeT(2),...
        pHingeT(3),Wg.Dihedral,[1;0;0],Wg.dRefPt');
end
% Rotate Panels with incidence angle
if Wg.Incidence ~= 0
    [Xw,Yw,Zw] = AxisRotating(Xw,Yw,Zw,Wg.Incidence,[0;1;0],Wg.iRefPt');
    [Xp,Yp,Zp] = AxisRotating(Xp,Yp,Zp,Wg.Incidence,[0;1;0],Wg.iRefPt');
    [pHingeR(1),pHingeR(2),pHingeR(3)] = AxisRotating(pHingeR(1),pHingeR(2),...
        pHingeR(3),Wg.Incidence,[0;1;0],Wg.iRefPt');
    [pHingeT(1),pHingeT(2),pHingeT(3)] = AxisRotating(pHingeT(1),pHingeT(2),...
        pHingeT(3),Wg.Incidence,[0;1;0],Wg.iRefPt');
end
pHinge = pHingeR';
rotAxis = pHingeT'-pHingeR';