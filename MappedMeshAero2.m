function [PanelDat,surfNum2]=MappedMeshAero2(AC,isurf,surfNum,PNP)
% Qinf=AC(1).State.Qinf;
% dt=AC(1).State.dt;
if AC(isurf).Config.Dihedral_ref == 1
    Pref = [0 0 0];% Rotate around x-axis
elseif AC(isurf).Config.Dihedral_ref == 2
    Pref = AC(isurf).Config.Proot;% master=man wing, slave=control surfaces
else
    error('please input ref point for dihedral rotation');
end
% lifting surface
Proot = AC(isurf).Config.Proot;%m
Ptip = AC(isurf).Config.Ptip;%m
Dihedral = AC(isurf).Config.Dihedral;
if Proot(2)==Ptip(2)
    surf_span = Ptip(:,3)-Proot(:,3);
    Ptip(:,3) = Proot(:,3);%m
    Ptip(:,2) = Ptip(:,2)+surf_span;%m
    Dihedral = Dihedral+90*pi/180;
end
rChord=AC(isurf).Config.rChord;%m
tChord=AC(isurf).Config.tChord;%m
X01=Proot(1);Y01=Proot(2);Z01=Proot(3);
X02=Ptip(1);Y02=Ptip(2);Z02=Ptip(3);
X03=X02+tChord;Y03=Y02;Z03=Z02;
X04=X01+rChord;Y04=Y01;Z04=Z01;

span=Ptip(2)-Proot(2);

% Grid generation
% nwake=1;% for steady flow
PanelDat.Nodes=[];
% PanelDat.NodesW=[];
PanelDat.VtxPt=[];
% PanelDat.Ring2Lift=[];
% PanelDat.VtxPtF=[];
% PanelDat.Vtx2panel=[];
PanelDat.ColPt=[];
PanelDat.NormV=[];
PanelDat.delx=[];
PanelDat.dely=[];
PanelDat.S=[];
PanelDat.WingPanel=[];
% PanelDat.WakePanel=[];
% PanelDat.TrailPanel=[];
PanelDat.SurfID=[];
PanelDat.SurfName=[];%Surface Name (Natee)
PanelDat.Panel_idx=[];%Panel index of surfaces (Natee)
PanelDat.rot_axis_vec=[];%Rotating axis in vector form for control surfaces (Natee)
PanelDat.rot_axis_pnt=[];%Rotating axis in point form for control surfaces (Natee)
%
nnode=[];npanel=[];%nwake=[];
surfNum2=0;
nspan_ass=[];
nchord_ass=[];
for i=1:length(AC(isurf).SubSurf)
    
    
    spanInt=AC(isurf).SubSurf(i).SpanData(1:2);
    nspan=AC(isurf).SubSurf(i).SpanData(3);
    chordInt=AC(isurf).SubSurf(i).ChordIntv;
    nchord=AC(isurf).SubSurf(i).ChordMesh;
    surfNum2=surfNum2+length(nchord);
    
    % identify sub-surface ID
    surfID=[];
    for j=1:length(nchord)
        surfID=[surfID;(surfNum+j)*ones(nchord(j),nspan)];
    end
%     surfNum=surfNum+surfNum2;
    surfNum=surfNum+length(nchord); %Fixed by Natee
    
    X1=X01+spanInt(1)*(X02-X01);
    X2=X01+spanInt(2)*(X02-X01);
    X3=X04+spanInt(2)*(X03-X04);
    X4=X04+spanInt(1)*(X03-X04);
    Y1=Y01+spanInt(1)*(Y02-Y01);
    Y2=Y01+spanInt(2)*(Y02-Y01);
    Y3=Y2;Y4=Y1;
    Z1=Z01;Z2=Z02;Z3=Z2;Z4=Z1;

    sj=linspace(chordInt(1),chordInt(2),nchord(1)+1);
    s0=sj;
    vt0=s0(1:end-1)+0.25*(s0(2:end)-s0(1:end-1));
    for j=2:length(nchord)
        sj=linspace(chordInt(j),chordInt(j+1),nchord(j)+1);
        s0=[s0 sj(2:end)];
        vt0=[vt0 sj(1:end-1)+0.25*(sj(2:end)-sj(1:end-1))];
    end
    vt0=[vt0 sj(end)+0.25*(sj(end)-sj(end-1))];
    %s0=-1+2*s0;% mapped to natural coordinates
    t0=linspace(-1,1,nspan+1);
    vt=-1+2*vt0;
    %[s1,t1]=meshgrid(s0,t0);s1=s1';t1=t1';% original wing grid
    [s2,t2]=meshgrid(vt,t0);s2=s2';t2=t2';% panel grid
    % Mapped from the natural coordinates to the real coordinates
    %x1=0.25*((1-s1).*(1-t1)*X1+(1-s1).*(1+t1)*X2+(1+s1).*(1+t1)*X3+(1+s1).*(1-t1)*X4);
    x2=0.25*((1-s2).*(1-t2)*X1+(1-s2).*(1+t2)*X2+(1+s2).*(1+t2)*X3+(1+s2).*(1-t2)*X4);
    
    %y1=0.25*((1-s1).*(1-t1)*Y1+(1-s1).*(1+t1)*Y2+(1+s1).*(1+t1)*Y3+(1+s1).*(1-t1)*Y4);
    y2=0.25*((1-s2).*(1-t2)*Y1+(1-s2).*(1+t2)*Y2+(1+s2).*(1+t2)*Y3+(1+s2).*(1-t2)*Y4);
    
    if strcmp(lower(AC(isurf).Config.Airfoil),'flat')
        z2=Proot(3)+zeros(size(x2));
    else % If there exists a wing camber
        Airfoil=load(AC(isurf).Config.Airfoil);
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
        Camber0=0.5*(UpperZ+LowerZ);
        delX=UpperX(end)-UpperX(1);% chord from the airfoil data
        %old interp (Ignore slop of the Last Panel)
%         Camber=interp1(UpperX+vt0(1),Camber0,vt0(1:end-1),'pchip','extrap')';%(Natee - Fix Camber interpolation + 1/4 panel chordwise distance (vt0(1)))
%         Camber=[Camber;Camber(end)+(vt0(end)-vt0(end-1))*(Camber(end-1)-Camber(end))/(vt0(end-2)-vt0(end-1))];
        %new interp (Interpolate all Panels)
        Camber=interp1(UpperX+vt0(1),Camber0,vt0(1:end),'spline')';%(Natee - Fix Camber interpolation + 1/4 panel chordwise distance (vt0(1)))
        %
        z2=Proot(3)+repmat(Camber,1,nspan+1).*...
            repmat(linspace(rChord,tChord,nspan+1),sum(nchord)+1,1)/delX;
    end
    panelDat=meshgen(x2,y2,z2,nspan,sum(nchord),...
        AC(isurf).Config.Incidence, Dihedral,...
        Pref,AC(isurf).Config.Symmetry,span);
%     panelDat=meshgen(x2,y2,z2,nspan,sum(nchord),1,...
%         AC(isurf).Config.Incidence, AC(isurf).Config.Dihedral,...
%         Pref,AC(1).State.Qinf,AC(1).State.dt,AC(isurf).Config.Symmetry,span);
    nvtx(i)=size(panelDat.VtxPt,1);
    nnode(i)=size(panelDat.Nodes,1);
    npanel(i)=size(panelDat.WingPanel,1);
%     nwake(i)=size(panelDat.NodesW,1);
    
    
    if strcmp(AC(isurf).Config.Symmetry,'yes')
        %Add surface name and panel index of control surfaces (Natee)
        PNP2=PNP+size(PanelDat.WingPanel,1);%Previous Npanel 
        %Right-side Surface
        PanelDat.SurfName=[PanelDat.SurfName,{[AC(isurf).SubSurf(i).Name 'R']}];
        PanelDat.Panel_idx=[PanelDat.Panel_idx,{PNP2+1:PNP2+(npanel(i))/2}];%Store Panel index of each surface (Natee)
        PanelDat.SurfID=[PanelDat.SurfID;reshape(surfID,nspan*sum(nchord),1)];
        
        %Left-side Surface
        PanelDat.SurfName=[PanelDat.SurfName,{[AC(isurf).SubSurf(i).Name 'L']}];
        PanelDat.Panel_idx=[PanelDat.Panel_idx,{PNP2+1+(npanel(i)/2):PNP2+(npanel(i))}];%Store Panel index of each surface (Natee)
        surfID=surfID+1;
        surfNum=surfNum+length(nchord);
        surfNum2=surfNum2+1;
        PanelDat.SurfID=[PanelDat.SurfID;reshape(surfID,nspan*sum(nchord),1)];
        
        %Add rotating axis for control surfaces (Natee)
        switch AC(isurf).SubSurf(i).Name
            case 'Flap'
                n2pR=panelDat.Nodes(panelDat.WingPanel(1,[1 2]),:);
                if n2pR(1,2)>n2pR(2,2)%Rotation Axis of the Right-side Flap pointing inward
                   n2pR=n2pR([2 1],:); 
                end
                n2pL=n2pR([2 1],:);%swap node 1-2 (flip)
                n2pL(:,2)=-n2pL(:,2);%Rotation Axis of the Left-side Flap pointing outward (flip + inverse Y)
            case 'Aileron'
                n2pR=panelDat.Nodes(panelDat.WingPanel(1,[1 2]),:);
                if n2pR(1,2)<n2pR(2,2)%Rotation Axis of the Right-side wing pointing outward
                   n2pR=n2pR([2 1],:); 
                end
                n2pL=n2pR;
                n2pL(:,2)=-n2pL(:,2);%Rotation Axis of the Left-side wing pointing outward (inverse Y)
            case 'Elevator'
                n2pR=panelDat.Nodes(panelDat.WingPanel(1,[1 2]),:);
                if n2pR(1,2)>n2pR(2,2)%Rotation Axis of the Right-side elevator pointing inward
                   n2pR=n2pR([2 1],:); 
                end
                n2pL=n2pR([2 1],:);%swap node 1-2 (flip)
                n2pL(:,2)=-n2pL(:,2);%Rotation Axis of the Left-side elevator pointing outward (flip + inverse Y)
            otherwise
                n2pR=[];
        end
        if numel(n2pR)>0
            vecR=n2pR(2,:)-n2pR(1,:);
            vecL=n2pL(2,:)-n2pL(1,:);
            PanelDat.rot_axis_pnt=[PanelDat.rot_axis_pnt,{n2pR}];
            PanelDat.rot_axis_pnt=[PanelDat.rot_axis_pnt,{n2pL}];
            PanelDat.rot_axis_vec=[PanelDat.rot_axis_vec,{vecR}];
            PanelDat.rot_axis_vec=[PanelDat.rot_axis_vec,{vecL}];
        else
            PanelDat.rot_axis_pnt=[PanelDat.rot_axis_pnt,{[]},{[]}];
            PanelDat.rot_axis_vec=[PanelDat.rot_axis_vec,{[]},{[]}];
        end
        %
    else
        PNP2=PNP+size(PanelDat.WingPanel,1);%Previous Npanel 
        PanelDat.SurfName=[PanelDat.SurfName,{AC(isurf).SubSurf(i).Name}];
        PanelDat.Panel_idx=[PanelDat.Panel_idx,{PNP2+1:PNP2+npanel(i)}];%Store Panel index of each surface (Natee)
        PanelDat.SurfID=[PanelDat.SurfID;reshape(surfID,nspan*sum(nchord),1)];
        
        %Add rotating axis for control surfaces (Natee)
        switch AC(isurf).SubSurf(i).Name
            case 'Rudder'
                n2p=panelDat.Nodes(panelDat.WingPanel(1,[1 2]),:);
                if n2p(1,3)<n2p(2,3)%Rotation Axis of the rudder is pointing upward
                   n2p=n2p([2 1],:); 
                end
            otherwise
                n2p=[];
        end
        if numel(n2p)>0
            rvec=n2p(2,:)-n2p(1,:);
            PanelDat.rot_axis_pnt=[PanelDat.rot_axis_pnt,{n2p}];
            PanelDat.rot_axis_vec=[PanelDat.rot_axis_vec,{rvec}];
        else
            PanelDat.rot_axis_pnt=[PanelDat.rot_axis_pnt,{[]}];
            PanelDat.rot_axis_vec=[PanelDat.rot_axis_vec,{[]}];
        end
        %
    end
    PanelDat.Nodes=[PanelDat.Nodes;panelDat.Nodes];
%     PanelDat.NodesW=[PanelDat.NodesW;panelDat.NodesW];
    PanelDat.VtxPt=[PanelDat.VtxPt;panelDat.VtxPt];
%     PanelDat.Ring2Lift=[PanelDat.Ring2Lift zeros(sum(npanel(1:i-1)),npanel(i));zeros(npanel(i),sum(npanel(1:i-1))) panelDat.Ring2Lift];
%     PanelDat.Vtx2panel=[PanelDat.Vtx2panel;panelDat.Vtx2panel+sum(nvtx(1:i-1))];
%     PanelDat.VtxPtF=[PanelDat.VtxPtF;panelDat.VtxPtF];
    PanelDat.ColPt=[PanelDat.ColPt;panelDat.ColPt];
    PanelDat.NormV=[PanelDat.NormV;panelDat.NormV];
    PanelDat.delx=[PanelDat.delx;panelDat.delx];
    PanelDat.dely=[PanelDat.dely;panelDat.dely];
    PanelDat.S=[PanelDat.S;panelDat.S];
    PanelDat.WingPanel=[PanelDat.WingPanel;panelDat.WingPanel+sum(nnode(1:i-1))];
%     PanelDat.WakePanel=[PanelDat.WakePanel;panelDat.WakePanel+sum(nwake(1:i-1))];
%     PanelDat.TrailPanel=[PanelDat.TrailPanel;panelDat.TrailPanel+sum(npanel(1:i-1))];
%     PanelDat.SurfID=[PanelDat.SurfID;reshape(surfID,nspan*sum(nchord),1)];
    
    if strcmp(AC(isurf).Config.Symmetry,'yes')
%         PanelDat.SurfID=[PanelDat.SurfID;reshape(surfID,nspan*sum(nchord),1)];
        nspan_ass=[nspan_ass;sum(nspan);sum(nspan)];%nspan of all subsurface
        nchord_ass=[nchord_ass;sum(nchord);sum(nchord)];%nchord of all subsurface
    else
        nspan_ass=[nspan_ass;sum(nspan)];%nspan of all subsurface
        nchord_ass=[nchord_ass;sum(nchord)];%nchord of all subsurface
    end
%     figure(5),clf,hold on
%     surf(x1,y1,0*x1,'facecolor','g','facealpha',0.5)
%     surf(x2,y2,z2,'facecolor','r','facealpha',0.5)
%     plot3([X1 X2 X3 X4 X1],[Y1 Y2 Y3 Y4 Y1],[0 0 0 0 0],'r-o')
%     xlabel('x'),ylabel('y')
%     pause
end
PanelDat.nvtx=sum(nvtx);
PanelDat.nnode=sum(nnode);
PanelDat.npanel=sum(npanel);
% PanelDat.nwake=sum(nwake);
PanelDat.nspan=nspan_ass;
PanelDat.nchord=nchord_ass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function panelDat=meshgen(x2,y2,z2,nspan,nchord,nwake,Incidence,Dihedral,...
%     Pref,Qinf,dt,Symmetry,span)
function panelDat=meshgen(x2,y2,z2,nspan,nchord,Incidence,Dihedral,...
    Pref,Symmetry,span)
nnode=(nchord+1)*(nspan+1);

% panel corner points
xx1=x2(1:end-1,1:end-1);
xx2=x2(1:end-1,2:end);
xx3=x2(2:end,2:end);
xx4=x2(2:end,1:end-1);

yy1=y2(1:end-1,1:end-1);
yy2=y2(1:end-1,2:end);
yy3=y2(2:end,2:end);
yy4=y2(2:end,1:end-1);

zz1=z2(1:end-1,1:end-1);
zz2=z2(1:end-1,2:end);
zz3=z2(2:end,2:end);
zz4=z2(2:end,1:end-1);

% panel areas
npanel=nchord*nspan;
dely=reshape(abs(yy2-yy1),npanel,1);
delx=reshape(abs(0.5*(xx3+xx4)-0.5*(xx1+xx2)),npanel,1);
S=reshape(dely.*delx,npanel,1);  

% % rotating due to Incidence
% if Incidence ~= 0
%     Xrot=Pref(1);Zrot=Pref(3);    
%     x22=Xrot+(x2-Xrot)*cos(Incidence)+(z2-Zrot)*sin(Incidence);
%     z2=Zrot-(x2-Xrot)*sin(Incidence)+(z2-Zrot)*cos(Incidence);
%     x2=x22;
% end
% 
% % rotating due to Dihedral
% if Dihedral ~= 0
%     Yrot=0;Zrot=0;% fuselage axis
%     y22=Yrot+(y2-Yrot)*cos(Dihedral)-(z2-Zrot)*sin(Dihedral);
%     z2=Zrot+(y2-Yrot)*sin(Dihedral)+(z2-Zrot)*cos(Dihedral);
%     y2=y22;
% end

%
% (Natee) Swap Dihedral with Incidence to avoid wing splitting issue
% rotating due to Dihedral
if Dihedral ~= 0
%     Yrot=0;Zrot=0;% fuselage axis
    Yrot=Pref(2);Zrot=Pref(3);%Natee fixed for winglet and offset Vtail gen
    %
    y22=Yrot+(y2-Yrot)*cos(Dihedral)-(z2-Zrot)*sin(Dihedral);
    z2=Zrot+(y2-Yrot)*sin(Dihedral)+(z2-Zrot)*cos(Dihedral);
    y2=y22;
end

% rotating due to Incidence
if Incidence ~= 0
    Xrot=Pref(1);Zrot=Pref(3);    
    x22=Xrot+(x2-Xrot)*cos(Incidence)+(z2-Zrot)*sin(Incidence);
    z2=Zrot-(x2-Xrot)*sin(Incidence)+(z2-Zrot)*cos(Incidence);
    x2=x22;
end
% 

%

% % wake vortices
% Xw1=x2(nchord+1,1);Yw1=y2(nchord+1,1);Zw1=z2(nchord+1,1);
% Xw2=x2(nchord+1,nspan+1);Yw2=y2(nchord+1,nspan+1);Zw2=z2(nchord+1,nspan+1);
% % Xw3=Xw2+Qinfx*dt;Yw3=Yw2;Zw3=Zw2;
% % Xw4=Xw1+Qinfx*dt;Yw4=Yw1;Zw4=Zw1;
% 
% % Natee -> wake = 6*span (Tornado's Default for validation)
% Lwake=6*span;
% wake_direction=Qinf/norm(Qinf);
% Xw3=Xw2+Lwake*wake_direction(1);
% Yw3=Yw2+Lwake*wake_direction(2);
% Zw3=Zw2+Lwake*wake_direction(3);
% 
% Xw4=Xw1+Lwake*wake_direction(1);
% Yw4=Yw1+Lwake*wake_direction(2);
% Zw4=Zw1+Lwake*wake_direction(3);
% 
% [sw,tw]=meshgrid(linspace(-1,1,nwake+1),linspace(-1,1,nspan+1));sw=sw';tw=tw';% wake
% xw=0.25*((1-sw).*(1-tw)*Xw1+(1-sw).*(1+tw)*Xw2+(1+sw).*(1+tw)*Xw3+(1+sw).*(1-tw)*Xw4);
% yw=0.25*((1-sw).*(1-tw)*Yw1+(1-sw).*(1+tw)*Yw2+(1+sw).*(1+tw)*Yw3+(1+sw).*(1-tw)*Yw4);
% zw=0.25*((1-sw).*(1-tw)*Zw1+(1-sw).*(1+tw)*Zw2+(1+sw).*(1+tw)*Zw3+(1+sw).*(1-tw)*Zw4);
%

% % wake panels
% nnum=cumsum(ones(nwake+1,nspan+1))+...
%     repmat((nwake+1)*(0:nspan),(nwake+1),1);
% nnum1=nnum(1:end-1,1:end-1);
% nnum2=nnum(1:end-1,2:end);
% nnum3=nnum(2:end,2:end);
% nnum4=nnum(2:end,1:end-1);
% WakeElem=[reshape(nnum1,nwake*nspan,1) reshape(nnum2,nwake*nspan,1) ...
%     reshape(nnum3,nwake*nspan,1) reshape(nnum4,nwake*nspan,1)];

% % corresponding trailing panels
% TrailElem=(nchord:nchord:npanel)';

% panel grid generation
npanel=nchord*nspan;
nnum=cumsum(ones(nchord+1,nspan+1))+...
    repmat((nchord+1)*(0:nspan),(nchord+1),1);
nnum1=nnum(1:end-1,1:end-1);
nnum2=nnum(1:end-1,2:end);
nnum3=nnum(2:end,2:end);
nnum4=nnum(2:end,1:end-1);
WingElem=[reshape(nnum1,npanel,1) reshape(nnum2,npanel,1) ...
    reshape(nnum3,npanel,1) reshape(nnum4,npanel,1)];

% panel corner points
xx1=x2(1:end-1,1:end-1);
xx2=x2(1:end-1,2:end);
xx3=x2(2:end,2:end);
xx4=x2(2:end,1:end-1);

yy1=y2(1:end-1,1:end-1);
yy2=y2(1:end-1,2:end);
yy3=y2(2:end,2:end);
yy4=y2(2:end,1:end-1);

zz1=z2(1:end-1,1:end-1);
zz2=z2(1:end-1,2:end);
zz3=z2(2:end,2:end);
zz4=z2(2:end,1:end-1);

% Vortex points
xv=reshape(0.5*(xx1+xx2),npanel,1);
yv=reshape(0.5*(yy1+yy2),npanel,1);
zv=reshape(0.5*(zz1+zz2),npanel,1);
% 
% enum=cumsum(ones(nchord,nspan))+...
%     repmat((nchord)*(0:nspan-1),(nchord),1);
% Ring2Lift=eye(npanel);
% iRing1=reshape(enum(2:end,:),(nchord-1)*nspan,1);
% iRing2=reshape(enum(1:end-1,:),(nchord-1)*nspan,1);
% Ring2Lift(iRing1,iRing2)=Ring2Lift(iRing1,iRing2)-eye(npanel-nspan);

% Vtx2panel
% nnode=(nchord+1)*(nspan+1);
% [nc1,ns1]=meshgrid(1:nspan,1:(nchord+1));
% numv1=(nchord+1)*(nc1-1)+ns1;
% 
% [nc2,ns2]=meshgrid(1:nspan+1,1:nchord);
% numv2=(nnode-(nchord+1))+nchord*(nc2-1)+ns2;
% Vtx2panel=[
%     reshape(numv1(1:end-1,:),npanel,1) ...
%     reshape(numv2(:,2:end),npanel,1) ...
%     reshape(numv1(2:end,:),npanel,1) ...
%     reshape(numv2(:,1:end-1),npanel,1)  
%     ];

% %Front Vortex points
% xvF=reshape(0.5*(xx1+xx2),npanel,1);
% yvF=reshape(0.5*(yy1+yy2),npanel,1);
% zvF=reshape(0.5*(zz1+zz2),npanel,1);

% Collocation points
xc=reshape(0.25*(xx1+xx2+xx3+xx4),npanel,1);
yc=reshape(0.25*(yy1+yy2+yy3+yy4),npanel,1);
zc=reshape(0.25*(zz1+zz2+zz3+zz4),npanel,1);

% normal vectors
vec1(:,:,1)=xx3-xx1;
vec1(:,:,2)=yy3-yy1;
vec1(:,:,3)=zz3-zz1;
vec2(:,:,1)=xx2-xx4;
vec2(:,:,2)=yy2-yy4;
vec2(:,:,3)=zz2-zz4;

nv=cross(vec1,vec2,3);
magnv=sqrt(dot(nv,nv,3));
nvx=nv(:,:,1)./magnv;
nvy=nv(:,:,2)./magnv;
nvz=nv(:,:,3)./magnv;
NormV=[reshape(nvx,npanel,1) reshape(nvy,npanel,1) reshape(nvz,npanel,1)];

% output
panelDat.Nodes=[reshape(x2,nnode,1) reshape(y2,nnode,1) reshape(z2,nnode,1)];
% panelDat.NodesW=[reshape(xw,2*(nspan+1),1) reshape(yw,2*(nspan+1),1) ...
%     reshape(zw,2*(nspan+1),1)];
panelDat.VtxPt=[xv yv zv];
% panelDat.Vtx2panel=Vtx2panel;
% panelDat.VtxPtF=[xvF yvF zvF];%
% panelDat.Ring2Lift=Ring2Lift;
panelDat.ColPt=[xc yc zc];
panelDat.NormV=NormV;
panelDat.delx=delx;
panelDat.dely=dely;
panelDat.S=S;
panelDat.WingPanel=WingElem;
% panelDat.WakePanel=WakeElem;
% panelDat.TrailPanel=TrailElem;

if strcmp(Symmetry,'yes')
    sidx1 = reshape(fliplr(reshape(nnode+(1:nnode),nchord+1,nspan+1)),[],1);
    sidx2 = reshape(fliplr(reshape(npanel+(1:npanel),nchord,nspan)),[],1);
    
    panelDat.Nodes=[panelDat.Nodes
        panelDat.Nodes(:,1) -panelDat.Nodes(:,2) panelDat.Nodes(:,3)
        ];
    panelDat.Nodes((nnode+1):end,:)=panelDat.Nodes(sidx1,:);
    

%     hspan=max(y2(1,:))-min(y2(:,1));
%     YW=panelDat.NodesW(:,2);
%     YW(1:2:end)=-YW(1:2:end);
%     YW(2:2:end)=flipud(YW(2:2:end)-hspan);
%     panelDat.NodesW=[panelDat.NodesW
%                      panelDat.NodesW(:,1) YW panelDat.NodesW(:,3)];
%     panelDat.NodesW(npanel+1:end,:)=panelDat.NodesW(sidx,:);
    
%     panelDat.NodesW=[panelDat.NodesW
%         panelDat.NodesW(:,1) -panelDat.NodesW(:,2) panelDat.NodesW(:,3)
%         ];

    %
    panelDat.VtxPt=[panelDat.VtxPt
        panelDat.VtxPt(:,1) -panelDat.VtxPt(:,2) panelDat.VtxPt(:,3)
        ];
    panelDat.VtxPt(npanel+1:end,:)=panelDat.VtxPt(sidx2,:);
    
%     panelDat.Ring2Lift=[Ring2Lift zeros(npanel)
%                         zeros(npanel) Ring2Lift];
%     panelDat.Vtx2panel=[Vtx2panel;nspan*(nchord+1)+nchord*(nspan+1)+Vtx2panel];
%     panelDat.VtxPtF=[panelDat.VtxPtF
%         panelDat.VtxPtF(:,1) -panelDat.VtxPtF(:,2) panelDat.VtxPtF(:,3)
%         ];
    panelDat.ColPt=[panelDat.ColPt
        panelDat.ColPt(:,1) -panelDat.ColPt(:,2) panelDat.ColPt(:,3)
        ];
    panelDat.ColPt(npanel+1:end,:)=panelDat.ColPt(sidx2,:);
    
    panelDat.NormV=[panelDat.NormV
        panelDat.NormV(:,1) -panelDat.NormV(:,2) panelDat.NormV(:,3)
        ];
    panelDat.NormV(npanel+1:end,:)=panelDat.NormV(sidx2,:);
    
    panelDat.delx=[panelDat.delx;panelDat.delx];
    panelDat.dely=[panelDat.dely;panelDat.dely];
    panelDat.S=[panelDat.S;panelDat.S];
    panelDat.WingPanel=[WingElem;nnode+WingElem];
%     panelDat.WakePanel=[WakeElem;(nwake+1)*(nspan+1)+WakeElem];
%     panelDat.TrailPanel=[TrailElem;npanel+TrailElem];
end

