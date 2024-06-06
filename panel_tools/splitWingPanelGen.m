function WgOut=splitWingPanelGen(Wg,iSurf)
%%Natee Add vertical Ptip input support
if Wg.Proot(2)==Wg.Ptip(2)
    Wg.Ptip(2)=Wg.Proot(2)+(Wg.Ptip(3)-Wg.Proot(3));
    Wg.Ptip(3)=Wg.Proot(3);
    Wg.Dihedral=90*pi/180;
    Wg.dRefPt=Wg.Proot;
end
%%

% Panel generation
if strcmp(Wg.Airfoil,'flat')
    Airfoil = [linspace(0,1,20)',zeros(20,1)];
else
    Airfoil=load(Wg.Airfoil);
end
if Airfoil(1,1)>0
    UpperX=Airfoil(2:Airfoil(1,1)+1,1);
    UpperZ=Airfoil(2:Airfoil(1,1)+1,2);
    LowerZ=Airfoil(Airfoil(1,1)+2:end,2);
else
    nAfoil=size(Airfoil,1);
    UpperX=Airfoil(1:nAfoil/2,1);
    UpperZ=Airfoil(1:nAfoil/2,2);
    LowerZ=Airfoil(nAfoil/2+1:end,2);
end

chordRatio=Wg.SubSurf(iSurf).ChordIntv(2);
nChordM=max(1,round(chordRatio*Wg.nChordPanel));%M for main wing
spStart=Wg.SubSurf(iSurf).SpanData(1);
spEnd=Wg.SubSurf(iSurf).SpanData(2);

[UpperX,nsort]=sort(UpperX);
UpperZ=UpperZ(nsort);
LowerZ=LowerZ(nsort);
CamberZ=0.5*(UpperZ+LowerZ);

% m for main wing, s for sub-wing
maxX=max(UpperX);minX=min(UpperX);
dX=maxX-minX;
delXm=0.25*(maxX-minX)/Wg.nChordPanel;
delXs=0.25*(1-chordRatio)*(maxX-minX)/Wg.SubSurf(iSurf).ChordMesh;
xmWing=linspace(minX,minX+chordRatio*dX,nChordM+1);
xsWing=linspace(minX+chordRatio*dX,maxX,Wg.SubSurf(iSurf).ChordMesh+1);
xWing=[xmWing xsWing(2:end)];
zWing=interp1(UpperX,CamberZ,xWing,'linear','extrap');
xmPanel=linspace(minX+delXm,(minX+chordRatio*dX)+delXm,nChordM+1);
xsPanel=linspace((minX+chordRatio*dX)+delXs,maxX+delXs,Wg.SubSurf(iSurf).ChordMesh+1);
xPanel=[xmPanel(1:end-1) xsPanel(1:end-1) maxX+delXm];
zPanel=interp1(UpperX,CamberZ,xPanel,'linear','extrap');
% figure(1),clf,hold on
% plot(UpperX,UpperZ,'-oc')
% plot(UpperX,LowerZ,'--')
% plot(xWing,zWing,'-sb')
% plot(xPanel,zPanel,':xr')
% pause

Wgm.nChordM=nChordM;
Wgm.xWing=(xWing-min(xWing))/(max(xWing)-min(xWing));
Wgm.zWing=zWing/(max(xWing)-min(xWing));
Wgm.xPanel=(xPanel-min(xPanel))/(max(xPanel)-min(xPanel));
Wgm.zPanel=zPanel/(max(xPanel)-min(xPanel));
Wgm.Proot=[Wg.Proot(1)+spStart*(Wg.Ptip(1)-Wg.Proot(1)),...
    Wg.Proot(2)+spStart*(Wg.Ptip(2)-Wg.Proot(2)),...
    Wg.Proot(3)+spStart*(Wg.Ptip(3)-Wg.Proot(3))];
Wgm.Ptip=[Wg.Proot(1)+spEnd*(Wg.Ptip(1)-Wg.Proot(1)),...
    Wg.Proot(2)+spEnd*(Wg.Ptip(2)-Wg.Proot(2)),...
    Wg.Proot(3)+spEnd*(Wg.Ptip(3)-Wg.Proot(3))];
Wgm.rChord=(Wg.rChord+spStart*(Wg.tChord-Wg.rChord));
Wgm.tChord=(Wg.rChord+spEnd*(Wg.tChord-Wg.rChord));
Wgm.chordRatio=chordRatio;
Wgm.dXroot=0.25*Wgm.rChord/Wg.nChordPanel;% leading root vortex ring position
Wgm.dXtip=0.25*Wgm.tChord/Wg.nChordPanel;
Wgm.nChordPanel=Wg.SubSurf(iSurf).ChordMesh+nChordM;
Wgm.nSpanPanel=max(1,round(abs(spEnd-spStart)*Wg.nSpanPanel));
Wgm.Incidence=Wg.Incidence;
Wgm.iRefPt=Wg.iRefPt;
Wgm.Dihedral=Wg.Dihedral;
Wgm.dRefPt=Wg.dRefPt;
[Xw,Yw,Zw,Xp,Yp,Zp,pHinge,rotAxis]=GenWingGrid(Wgm);%w for wings, p for panels

pnNum=reshape(1:Wgm.nChordPanel*Wgm.nSpanPanel,Wgm.nChordPanel,Wgm.nSpanPanel);
subSurfNum=reshape(pnNum(nChordM+1:Wgm.nChordPanel,:),1,...
    Wg.SubSurf(iSurf).ChordMesh*Wgm.nSpanPanel);
WgOut.RpanelNum=subSurfNum;%sub-surface panel numbers
WgOut.HingePt_R=pHinge;% sub-surface hinge position
WgOut.RotAxis_R=rotAxis;% sub-surface rotating axis

arrow3dRoundHead(pHinge',pHinge'+rotAxis')

WgOut.Right.Xw=Xw;WgOut.Right.Yw=Yw;WgOut.Right.Zw=Zw;
WgOut.Right.Xp=Xp;WgOut.Right.Yp=Yp;WgOut.Right.Zp=Zp;
if strcmp(Wg.Symmetry,'yes')
    WgOut.Left.Xw=Xw;WgOut.Left.Yw=-Yw;WgOut.Left.Zw=Zw;
    WgOut.Left.Xp=Xp;WgOut.Left.Yp=-Yp;WgOut.Left.Zp=Zp;

    WgOut.LpanelNum=Wgm.nChordPanel*Wgm.nSpanPanel+...
        subSurfNum;
    WgOut.HingePt_L=[pHinge(1) -pHinge(2) pHinge(3)];% sub-surface hinge position
    WgOut.RotAxis_L=[rotAxis(1) -rotAxis(2) rotAxis(3)];% sub-surface rotating axis
    arrow3dRoundHead(WgOut.HingePt_L',WgOut.HingePt_L'+WgOut.RotAxis_L')
end

if iSurf==1;sclr=0.7*ones(1,3);else;sclr=0.9*ones(1,3);end
mesh(Xw(1:nChordM+1,:),Yw(1:nChordM+1,:),Zw(1:nChordM+1,:),'facecolor','none','EdgeColor','r')
surf(Xp(1:nChordM+1,:),Yp(1:nChordM+1,:),Zp(1:nChordM+1,:),'facecolor','g','edgecolor','none','facealpha',0.75)
mesh(Xw(nChordM+1:end,:),Yw(nChordM+1:end,:),Zw(nChordM+1:end,:),'facecolor','none','EdgeColor','b')
surf(Xp(nChordM+1:end,:),Yp(nChordM+1:end,:),Zp(nChordM+1:end,:),'facecolor',sclr,'edgecolor','none','facealpha',0.75)

if strcmp(Wg.Symmetry,'yes')
    mesh(Xw(1:nChordM+1,:),-Yw(1:nChordM+1,:),Zw(1:nChordM+1,:),'facecolor','none','EdgeColor','r')
    surf(Xp(1:nChordM+1,:),-Yp(1:nChordM+1,:),Zp(1:nChordM+1,:),'facecolor','g','edgecolor','none','facealpha',0.75)
    mesh(Xw(nChordM+1:end,:),-Yw(nChordM+1:end,:),Zw(nChordM+1:end,:),'facecolor','none','EdgeColor','b')
    surf(Xp(nChordM+1:end,:),-Yp(nChordM+1:end,:),Zp(nChordM+1:end,:),'facecolor',sclr,'edgecolor','none','facealpha',0.75)
end
0;
%%%%%%%%%%% sub-functions %%%%%%%%%%%
function [Xw,Yw,Zw,Xp,Yp,Zp,pHinge,rotAxis]=GenWingGrid(Wg)
ProotW=Wg.Proot;
ProotP=[Wg.Proot(1)+Wg.dXroot Wg.Proot(2:3)];
PtipW=Wg.Ptip;
PtipP=[Wg.Ptip(1)+Wg.dXtip Wg.Ptip(2:3)];

% pHingeR=[ProotW(1)+Wg.chordRatio*Wg.rChord;ProotW(2);ProotW(3)];
% pHingeT=[PtipW(1)+Wg.chordRatio*Wg.tChord;PtipW(2);PtipW(3)];

Wg.nChordPanel=Wg.nChordPanel+1;% nChord intervals with nChord+1 nodes
Wg.nSpanPanel=Wg.nSpanPanel+1;% nSpan intervals with nSpan+1 nodes
x0W=linspace(ProotW(1),PtipW(1),Wg.nSpanPanel);
x0P=linspace(ProotP(1),PtipP(1),Wg.nSpanPanel);
yi=linspace(Wg.Proot(2),Wg.Ptip(2),Wg.nSpanPanel);
chDist=linspace(Wg.rChord,Wg.tChord,Wg.nSpanPanel);
for i=1:Wg.nChordPanel
    Yw(i,:)=yi;
    Xw(i,:)=x0W+Wg.xWing(i)*chDist;
    Zw(i,:)=ProotW(3)+Wg.zWing(i)*chDist;
    Yp(i,:)=Yw(i,:);
    Xp(i,:)=x0P+Wg.xPanel(i)*chDist;
    Zp(i,:)=ProotP(3)+Wg.zPanel(i)*chDist;
end

% Natee - Fix Hinge position -> move Hinge from structural grid to aerodynamic grid
pHingeR=[Xp(Wg.nChordM+1,1);Yp(Wg.nChordM+1,1);Zp(Wg.nChordM+1,1)];
pHingeT=[Xp(Wg.nChordM+1,end);Yp(Wg.nChordM+1,end);Zp(Wg.nChordM+1,end)];
%

%Edited swap dihedral and incidence rotation to preven trailing edge separation at root (Natee)

% rotating due to Dihedral
% if Wg.Dihedral ~= 0
%     Yw2=Wg.dRefPt(2)+(Yw-Wg.dRefPt(2))*cos(Wg.Dihedral)-(Zw-Wg.dRefPt(3))*sin(Wg.Dihedral);
%     Zw=Wg.dRefPt(3)+(Yw-Wg.dRefPt(2))*sin(Wg.Dihedral)+(Zw-Wg.dRefPt(3))*cos(Wg.Dihedral);
%     Yw=Yw2;
%     Yp2=Wg.dRefPt(2)+(Yp-Wg.dRefPt(2))*cos(Wg.Dihedral)-(Zp-Wg.dRefPt(3))*sin(Wg.Dihedral);
%     Zp=Wg.dRefPt(3)+(Yp-Wg.dRefPt(2))*sin(Wg.Dihedral)+(Zp-Wg.dRefPt(3))*cos(Wg.Dihedral);
%     Yp=Yp2;
% end
if Wg.Dihedral ~= 0
    [Xw,Yw,Zw] = AxisRotating(Xw,Yw,Zw,Wg.Dihedral,[1;0;0],Wg.dRefPt');
    [Xp,Yp,Zp] = AxisRotating(Xp,Yp,Zp,Wg.Dihedral,[1;0;0],Wg.dRefPt');
    [pHingeR(1),pHingeR(2),pHingeR(3)] = AxisRotating(pHingeR(1),pHingeR(2),...
        pHingeR(3),Wg.Dihedral,[1;0;0],Wg.dRefPt');
    [pHingeT(1),pHingeT(2),pHingeT(3)] = AxisRotating(pHingeT(1),pHingeT(2),...
        pHingeT(3),Wg.Dihedral,[1;0;0],Wg.dRefPt');
end

% rotating due to Incidence
% if Wg.Incidence ~= 0
%     Xw2=Wg.iRefPt(1)+(Xw-Wg.iRefPt(1))*cos(Wg.Incidence)+(Zw-Wg.iRefPt(3))*sin(Wg.Incidence);
%     Zw=Wg.iRefPt(3)-(Xw-Wg.iRefPt(1))*sin(Wg.Incidence)+(Zw-Wg.iRefPt(3))*cos(Wg.Incidence);
%     Xw=Xw2;
%     Xp2=Wg.iRefPt(1)+(Xp-Wg.iRefPt(1))*cos(Wg.Incidence)+(Zp-Wg.iRefPt(3))*sin(Wg.Incidence);
%     Zp=Wg.iRefPt(3)-(Xp-Wg.iRefPt(1))*sin(Wg.Incidence)+(Zp-Wg.iRefPt(3))*cos(Wg.Incidence);
%     Xp=Xp2;
% end
if Wg.Incidence ~= 0
    [Xw,Yw,Zw] = AxisRotating(Xw,Yw,Zw,Wg.Incidence,[0;1;0],Wg.iRefPt');
    [Xp,Yp,Zp] = AxisRotating(Xp,Yp,Zp,Wg.Incidence,[0;1;0],Wg.iRefPt');
    [pHingeR(1),pHingeR(2),pHingeR(3)] = AxisRotating(pHingeR(1),pHingeR(2),...
        pHingeR(3),Wg.Incidence,[0;1;0],Wg.iRefPt');
    [pHingeT(1),pHingeT(2),pHingeT(3)] = AxisRotating(pHingeT(1),pHingeT(2),...
        pHingeT(3),Wg.Incidence,[0;1;0],Wg.iRefPt');
end

pHinge=pHingeR';
rotAxis=pHingeT'-pHingeR';


