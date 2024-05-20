function [AC, PanelDat, AcControl, AcTrim]=PanelGen04(inputfile,state)
[AC, AcControl, AcTrim]=feval(inputfile,state);
nSurf = 0;
for i=1:numel(AC)
    try
        AC(i).SubSurf(1);
        nSurf = nSurf+1;
    end
end
PanelDat.Nodes=[];
PanelDat.NodesW=[];
PanelDat.VtxPt=[];%All Vortex points (4 points at middle of vortex lines of each ring)
% PanelDat.Ring2Lift=[];
% PanelDat.Vtx2panel=[];
% PanelDat.VtxPtF=[];%Vortex points at the front of each ring
PanelDat.ColPt=[];
PanelDat.NormV=[];
PanelDat.delx=[];
PanelDat.dely=[];
PanelDat.S=[];
PanelDat.WingPanel=[];
% PanelDat.WakePanel=[];
% PanelDat.TrailPanel=[];
PanelDat.SurfID=[];
nnode=[];npanel=[];%nwake=[];
PanelDat.nspan=[];
PanelDat.nchord=[];
surfNum1=0;
%
PanelDat.SurfName=[];%Surface Name (Natee)
PanelDat.Panel_idx=[];%Panel index of surfaces (Natee)
PanelDat.rot_axis_vec=[];%Rotating axis in vector form for control surfaces (Natee)
PanelDat.rot_axis_pnt=[];%Rotating axis in 2 points form for control surfaces (Natee)
%
NP=zeros(1,nSurf);

for i=1:nSurf
    PNP=size(PanelDat.WingPanel,1);%Previous Npanel (Natee)
    
    [panelDat,surfNum2]=MappedMeshAero2(AC,i,surfNum1,PNP);%Send Previous Npanel for surface indexing (Natee)
    surfNum1=surfNum1+surfNum2;
    
    nvtx(i)=size(panelDat.VtxPt,1);
    nnode(i)=size(panelDat.Nodes,1);
    npanel(i)=size(panelDat.WingPanel,1);
%     nwake(i)=size(panelDat.NodesW,1);

    PanelDat.Nodes=[PanelDat.Nodes;panelDat.Nodes];
%     PanelDat.NodesW=[PanelDat.NodesW;panelDat.NodesW];
    PanelDat.VtxPt=[PanelDat.VtxPt;panelDat.VtxPt];
    %Fix Ring2Lift for multiple Surf model
%     if sum(NP)==0
%         PanelDat.Ring2Lift=panelDat.Ring2Lift;
%     else
%         idx1=NP(i-1)+1;
%         idx2=NP(i-1)+size(panelDat.Ring2Lift,1);
%         PanelDat.Ring2Lift(idx1:idx2,idx1:idx2)=panelDat.Ring2Lift;
%     end
%     NP(i)=size(PanelDat.Ring2Lift,1);

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
    PanelDat.SurfID=[PanelDat.SurfID;panelDat.SurfID];
    
    PanelDat.nspan=[PanelDat.nspan;panelDat.nspan];
    PanelDat.nchord=[PanelDat.nchord;panelDat.nchord];
    
    PanelDat.SurfName=[PanelDat.SurfName,panelDat.SurfName];%Surface Name (Natee)
    PanelDat.Panel_idx=[PanelDat.Panel_idx,panelDat.Panel_idx];%Panel index of surfaces (Natee)
    PanelDat.rot_axis_vec=[PanelDat.rot_axis_vec,panelDat.rot_axis_vec];%Rotating axis in vector form for control surfaces (Natee)
    PanelDat.rot_axis_pnt=[PanelDat.rot_axis_pnt,panelDat.rot_axis_pnt];%Rotating axis in 2 points form for control surfaces (Natee)
end
PanelDat.nvtx=nvtx;
PanelDat.nnode=nnode;
PanelDat.npanel=npanel;
% PanelDat.nwake=nwake;

%Get Additional Panel Data for Lattice Generation

%Get number of panels of all main surfaces
count=0;
for i=1:numel(AC)
    npaneli = [];
    if numel(AC(i).SubSurf) > 0
        count = count+1;
        for j=1:numel(AC(i).SubSurf)
            npaneli(j) = AC(i).SubSurf(j).SpanData(3)*sum(AC(i).SubSurf(j).ChordMesh);
        end
        
        if strcmp(AC(i).Config.Symmetry,'yes');
            npanel_mainsurf(count) = 2*sum(npaneli);
        else
            npanel_mainsurf(count) = sum(npaneli);
        end
    end
end

%Get Trailing edge position
count=0;
for i=1:numel(AC)
    if numel(AC(i).SubSurf) > 0
        count = count+1;
        
%         pdx = (AC(i).SubSurf(1).ChordIntv(2)-AC(i).SubSurf(1).ChordIntv(1))/AC(i).SubSurf(1).ChordMesh;
%         dx_root = 0.25*pdx*AC(i).Config.rChord;
%         dx_tip = 0.25*pdx*AC(i).Config.tChord;
        
        TE = zeros(2,3);
        TE(1,:) = AC(i).Config.Proot;
        TE(2,:) = AC(i).Config.Ptip;
        TE(1,1) = TE(1,1)+AC(i).Config.rChord;%+dx_root;
        TE(2,1) = TE(2,1)+AC(i).Config.tChord;%+dx_tip;
        
        Pref = AC(i).Config.Proot;
        Incidence = AC(i).Config.Incidence;
        Dihedral = AC(i).Config.Dihedral;
        
        
        
        xx = TE(:,1);
        yy = TE(:,2);
        zz = TE(:,3);
        
        %rotating due to Dihedral
        if Dihedral ~= 0
%             Yrot=0;Zrot=0;% fuselage axis
            Yrot=Pref(2);Zrot=Pref(3);%Natee fixed for winglet and offset Vtail gen
            %
            yy2=Yrot+(yy-Yrot)*cos(Dihedral)-(zz-Zrot)*sin(Dihedral);
            zz2=Zrot+(yy-Yrot)*sin(Dihedral)+(zz-Zrot)*cos(Dihedral);
            xx2=xx;
            TE=[xx2,yy2,zz2];
            xx=xx2;
            yy=yy2;
            zz=zz2;
        end
        
        %rotating due to Incidence
        if Incidence ~= 0
            Xrot=Pref(1);Zrot=Pref(3);
            xx2=Xrot+(xx-Xrot)*cos(Incidence)+(zz-Zrot)*sin(Incidence);
            zz2=Zrot-(xx-Xrot)*sin(Incidence)+(zz-Zrot)*cos(Incidence);
            yy2=yy;
            TE=[xx2,yy2,zz2];
%             xx=xx2;
%             yy=yy2;
%             zz=zz2;
        end
        
        Trailing_Edge{count} = TE;
        TE_mirror = TE;
        TE_mirror(:,2) = -TE(:,2);
        Trailing_Edge_mirror{count} = TE_mirror;
    end
end

npanel_mainsurf2 = cumsum(npanel_mainsurf);
main_panel_idx = zeros(sum(npanel_mainsurf),1);
for i=1:numel(npanel_mainsurf2)
    if i==1
        main_panel_idx(1:npanel_mainsurf2(i)) = i;
    else
        main_panel_idx((npanel_mainsurf2(i-1)+1):npanel_mainsurf2(i)) = i;
    end
end
PanelDat.main_panel_idx = main_panel_idx;
PanelDat.Trailing_Edge = Trailing_Edge;
PanelDat.Trailing_Edge_mirror = Trailing_Edge_mirror;
% 
% figure(999);hold on;
% plot3(PanelDat.Nodes(:,1),PanelDat.Nodes(:,2),PanelDat.Nodes(:,3),'b.');
% for i=1:numel(PanelDat.Trailing_Edge)
% plot3(PanelDat.Trailing_Edge{i}(:,1),PanelDat.Trailing_Edge{i}(:,2),PanelDat.Trailing_Edge{i}(:,3),'r-');
% end
% axis equal;
0;