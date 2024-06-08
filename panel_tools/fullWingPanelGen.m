function WgOut=fullWingPanelGen(Wg,SpanData)
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
spStart=SpanData(1);
spEnd=SpanData(2);
[UpperX,nsort]=sort(UpperX);
UpperZ=UpperZ(nsort);
LowerZ=LowerZ(nsort);
CamberZ=0.5*(UpperZ+LowerZ);
delX0=0.25*(max(UpperX)-min(UpperX))/Wg.nChordPanel;
xWing=linspace(min(UpperX),max(UpperX),Wg.nChordPanel+1);
zWing=interp1(UpperX,CamberZ,xWing,'linear','extrap');
xPanel=linspace(min(UpperX)+delX0,max(UpperX)+delX0,Wg.nChordPanel+1);
zPanel=interp1(UpperX,CamberZ,xPanel,'linear','extrap');
% figure(1),clf,hold on
% plot(UpperX,UpperZ,'-oc')
% plot(UpperX,LowerZ,'--')
% plot(UpperX,CamberZ,'--')
% plot(xWing,zWing,'sb')
% plot(xPanel,zPanel,'xr')
% pause

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
Wgm.nChordPanel=Wg.nChordPanel;
Wgm.nSpanPanel=max(1,round(abs(spEnd-spStart)*Wg.nSpanPanel));
Wgm.Incidence=Wg.Incidence;
Wgm.iRefPt=Wg.iRefPt;
Wgm.Dihedral=Wg.Dihedral;
Wgm.dRefPt=Wg.dRefPt;
[Xw,Yw,Zw,Xp,Yp,Zp]=GenWingGrid(Wgm);%w for wings, p for panels

WgOut.Right.Xw=Xw;WgOut.Right.Yw=Yw;WgOut.Right.Zw=Zw;
WgOut.Right.Xp=Xp;WgOut.Right.Yp=Yp;WgOut.Right.Zp=Zp;
if strcmp(Wg.Symmetry,'yes')
    WgOut.Left.Xw=Xw;WgOut.Left.Yw=-Yw;WgOut.Left.Zw=Zw;
    WgOut.Left.Xp=Xp;WgOut.Left.Yp=-Yp;WgOut.Left.Zp=Zp;
end
mesh(Xw,Yw,Zw,'facecolor','none','EdgeColor','r')
surf(Xp,Yp,Zp,'facecolor','g','facealpha',0.75)
if strcmp(Wg.Symmetry,'yes')
    mesh(Xw,-Yw,Zw,'facecolor','none','EdgeColor','r')
    surf(Xp,-Yp,Zp,'facecolor','g','facealpha',0.75)
end
%%%%%%%%%%% sub-functions %%%%%%%%%%%
function [Xw,Yw,Zw,Xp,Yp,Zp]=GenWingGrid(Wg)
dXroot=0.25*Wg.rChord/Wg.nChordPanel;% leading root vortex ring position
dXtip=0.25*Wg.tChord/Wg.nChordPanel;
ProotW=Wg.Proot;
ProotP=[Wg.Proot(1)+dXroot Wg.Proot(2:3)];
PtipW=Wg.Ptip;
PtipP=[Wg.Ptip(1)+dXtip Wg.Ptip(2:3)];

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
end

