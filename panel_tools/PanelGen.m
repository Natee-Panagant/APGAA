function [PanelDat,FC,Sc,Sm,Si,So,S,pspan,pchord,normvec] = PanelGen(input_filename)
%%%%%%%%%%%%%%%%%%%%
% Panel Generation %
%%%%%%%%%%%%%%%%%%%%
% Read input file
[AC,FC] = feval(input_filename);
nSurf = size(AC,2);

PanelDat.Nodes = [];
PanelDat.NodesW = [];
PanelDat.VtxPt = [];
PanelDat.Ring2Lift = [];
PanelDat.ColPt = [];
PanelDat.NormV = [];
PanelDat.delx = [];
PanelDat.dely = [];
PanelDat.S = [];
PanelDat.WingPanel = [];
PanelDat.WakePanel = [];
PanelDat.TrailPanel = [];
PanelDat.SurfID = [];
nnode = [];nwake = [];npanel = [];

iSurf = 0;
iSS = 0;% Sub-Surface counter
for i=1:nSurf
    iLabel1 = strcmp({AC(i).Label}, 'Wing');%Wing
    iLabel2 = strcmp({AC(i).Label}, 'Horizontal tail');%Horizontal tail
    iLabel3 = strcmp({AC(i).Label}, 'Vertical tail');%Vertical tail
    iLabel4 = strcmp({AC(i).Label}, 'Fuselage');%Vertical tail

    if iLabel1||iLabel2||iLabel3||iLabel4
        % Add SubSurf field if not input
        if ~isfield(AC(i),'SubSurf')
            AC(i).SubSurf = [];
        end
        if size(AC(i).SubSurf,2) == 0 && numel(AC(i).nSpanPanel) > 0
            iSurf = iSurf+1;
            AcSurf = fullWingPanelGen(AC(i),[0 1]);
            panelDat = MappedMeshAero(AC(i),AcSurf,2);%second mode panel index see details in MappedMeshAero.m

            nvtx(iSurf) = size(panelDat.VtxPt,1);
            nnode(iSurf) = size(panelDat.Nodes,1);
            npanel(iSurf) = size(panelDat.WingPanel,1);
            nwake(iSurf) = size(panelDat.NodesW,1);

            %Natee - Store Name of Sub-Surfaces and their panel index
            iSS = iSS+1;
            NPi = size(panelDat.WingPanel,1);
            PanelDat.SurfName{iSS} = AC(i).Label; %Name of all sub-surfaces
            PanelDat.SurfPanelID{iSS} = sum(npanel(1:iSurf-1))+1:sum(npanel(1:iSurf-1))+NPi;
            PanelDat.HingePt{iSS} = [];%neglect for non-control surface
            PanelDat.HingeAxis{iSS} = [];%neglect for non-control surface
            PanelDat.RotAngle{iSS} = [];%neglect for non-control surface
            %

            PanelDat.Nodes = [PanelDat.Nodes;panelDat.Nodes];
            PanelDat.NodesW = [PanelDat.NodesW;panelDat.NodesW];
            PanelDat.VtxPt = [PanelDat.VtxPt;panelDat.VtxPt];
            PanelDat.Ring2Lift = [PanelDat.Ring2Lift zeros(sum(npanel(1:iSurf-1)),npanel(iSurf));zeros(npanel(iSurf),sum(npanel(1:iSurf-1))) panelDat.Ring2Lift];
            PanelDat.ColPt = [PanelDat.ColPt;panelDat.ColPt];
            PanelDat.NormV = [PanelDat.NormV;panelDat.NormV];
            PanelDat.delx = [PanelDat.delx;panelDat.delx];
            PanelDat.dely = [PanelDat.dely;panelDat.dely];
            PanelDat.S = [PanelDat.S;panelDat.S];
            PanelDat.WingPanel = [PanelDat.WingPanel;panelDat.WingPanel+sum(nnode(1:iSurf-1))];
            PanelDat.WakePanel = [PanelDat.WakePanel;panelDat.WakePanel+sum(nwake(1:iSurf-1))];
            PanelDat.TrailPanel = [PanelDat.TrailPanel;panelDat.TrailPanel+sum(npanel(1:iSurf-1))];
        elseif size(AC(i).SubSurf,2)>0
            [iInterval,iSubSurf] = getWingSec(AC(i));
            NP0 = sum(npanel);
            NPi = 0;
            init = 0;
            k = 0;
            for j = 1:length(iSubSurf)
                iSurf = iSurf+1;
                if iSubSurf(j) == 0
                    AcSurf = fullWingPanelGen(AC(i),[iInterval(j) iInterval(j+1)]);
                    panelDat = MappedMeshAero(AC(i),AcSurf,2);%second mode panel index see details in MappedMeshAero.m
                    if init == 0
                        k = k+1;
                        SSName{k} = AC(i).Label;
                        SSPidx{k} = [];%add later after loop complete
                        SSHingePt{k} = [];%neglect for non-control surface
                        SSHingeAxis{k} = [];%neglect for non-control surface
                        SSRotAngle{k} = [];%neglect for non-control surface
                        init = 1;
                    end
                    NPi = NPi+size(panelDat.WingPanel,1);
                else
                    AcSurf = splitWingPanelGen(AC(i),iSubSurf(j));
                    iPanel = sum(npanel(1:iSurf-1));
                    AC(i).SubSurf(iSubSurf(j)).RpanelNum = iPanel+AcSurf.RpanelNum;
                    AC(i).SubSurf(iSubSurf(j)).HingePt_R = AcSurf.HingePt_R;
                    AC(i).SubSurf(iSubSurf(j)).RotAxis_R = AcSurf.RotAxis_R;
                    if strcmp(AC(i).Symmetry,'yes')
                        AC(i).SubSurf(iSubSurf(j)).LpanelNum = iPanel+AcSurf.LpanelNum;
                        AC(i).SubSurf(iSubSurf(j)).HingePt_L = AcSurf.HingePt_L;
                        AC(i).SubSurf(iSubSurf(j)).RotAxis_L = AcSurf.RotAxis_L;
                    end
                    panelDat = MappedMeshAero(AC(i),AcSurf,2);%second mode panel index see details in MappedMeshAero.m
                    if numel(AcSurf.RpanelNum)<size(panelDat.WingPanel,1) && init == 0
                        k = k+1;
                        SSName{k} = AC(i).Label;
                        SSPidx{k} = [];%add later after loop complete
                        SSHingePt{k} = [];%neglect for non-control surface
                        SSHingeAxis{k} = [];%neglect for non-control surface
                        SSRotAngle{k} = [];%neglect for non-control surface
                        init = 1;
                    end

                    if strcmp(AC(i).Symmetry,'yes')
                        k = k+1;
                        SSName{k} = [AC(i).SubSurf(iSubSurf(j)).Name '_R'];
                        SSPidx{k} = iPanel+AcSurf.RpanelNum;
                        SSHingePt{k} = AcSurf.HingePt_R;%Hinge Point
                        SSHingeAxis{k} = AcSurf.RotAxis_R;%Hinge Axis
                        SSRotAngle{k} = AC(i).SubSurf(iSubSurf(j)).RotAngle;%Rotation angle of control surfaces

                        k = k+1;
                        SSName{k} = [AC(i).SubSurf(iSubSurf(j)).Name '_L'];
                        SSPidx{k} = iPanel+AcSurf.LpanelNum;
                        SSHingePt{k} = AcSurf.HingePt_L;
                        SSHingeAxis{k} = AcSurf.RotAxis_L;
                        SSRotAngle{k} = AC(i).SubSurf(iSubSurf(j)).RotAngle;%Rotation angle of control surfaces
                    else
                        k=k+1;
                        SSName{k} = AC(i).SubSurf(iSubSurf(j)).Name;
                        SSPidx{k} = iPanel+AcSurf.RpanelNum;
                        SSHingePt{k} = AcSurf.HingePt_R;%Hinge Point
                        SSHingeAxis{k} = AcSurf.RotAxis_R;%Hinge Axis
                        SSRotAngle{k} = AC(i).SubSurf(iSubSurf(j)).RotAngle;%Rotation angle of control surfaces
                    end
                    NPi = NPi+size(panelDat.WingPanel,1);
                end

                nvtx(iSurf) = size(panelDat.VtxPt,1);
                nnode(iSurf) = size(panelDat.Nodes,1);
                npanel(iSurf) = size(panelDat.WingPanel,1);
                nwake(iSurf) = size(panelDat.NodesW,1);
                PanelDat.Nodes = [PanelDat.Nodes;panelDat.Nodes];
                PanelDat.NodesW = [PanelDat.NodesW;panelDat.NodesW];
                PanelDat.VtxPt = [PanelDat.VtxPt;panelDat.VtxPt];
                PanelDat.Ring2Lift = [PanelDat.Ring2Lift zeros(sum(npanel(1:iSurf-1)),npanel(iSurf));zeros(npanel(iSurf),sum(npanel(1:iSurf-1))) panelDat.Ring2Lift];
                PanelDat.ColPt = [PanelDat.ColPt;panelDat.ColPt];
                PanelDat.NormV = [PanelDat.NormV;panelDat.NormV];
                PanelDat.delx = [PanelDat.delx;panelDat.delx];
                PanelDat.dely = [PanelDat.dely;panelDat.dely];
                PanelDat.S = [PanelDat.S;panelDat.S];
                PanelDat.WingPanel = [PanelDat.WingPanel;panelDat.WingPanel+sum(nnode(1:iSurf-1))];
                PanelDat.WakePanel = [PanelDat.WakePanel;panelDat.WakePanel+sum(nwake(1:iSurf-1))];
                PanelDat.TrailPanel = [PanelDat.TrailPanel;panelDat.TrailPanel+sum(npanel(1:iSurf-1))];
            end
            %Natee - Store Name of Sub-Surfaces and their panel index
            init = 0;
            for k = 1:numel(SSName)
                if numel(SSPidx{k})>0
                    iSS = iSS+1;
                    PanelDat.SurfName{iSS} = SSName{k}; %Name of all sub-surfaces
                    PanelDat.SurfPanelID{iSS} = SSPidx{k};
                    PanelDat.HingePt{iSS} = SSHingePt{k};
                    PanelDat.HingeAxis{iSS} = SSHingeAxis{k};
                    PanelDat.RotAngle{iSS} = SSRotAngle{k};
                else
                    if init == 0
                        iSS = iSS+1;
                        PanelDat.SurfName{iSS} = SSName{k}; %Name of all sub-surfaces
                        PanelDat.SurfPanelID{iSS} = setdiff(NP0+(1:NPi),cell2mat(SSPidx));
                        PanelDat.HingePt{iSS} = SSHingePt{k};
                        PanelDat.HingeAxis{iSS} = SSHingeAxis{k};
                        PanelDat.RotAngle{iSS} = SSRotAngle{k};
                        init = 1;
                    end
                end
            end
            clear SSName SSPidx SSHingePt SSHingeAxis SSRotAngle
            %
        end
    end
end
PanelDat.nvtx = nvtx;
PanelDat.nnode = nnode;
PanelDat.npanel = npanel;
PanelDat.nwake = nwake;


Trl = PanelDat.Ring2Lift;
nPanel = size(PanelDat.WingPanel,1);

% Rotate Control-Surfaces %Natee
for i = 1:numel(PanelDat.RotAngle)
    if PanelDat.RotAngle{i} ~= 0
        eidx = PanelDat.SurfPanelID{i};
        nidx = unique(reshape(PanelDat.WingPanel(eidx,:),[],1));
        X0 = PanelDat.Nodes(nidx,1);
        Y0 = PanelDat.Nodes(nidx,2);
        Z0 = PanelDat.Nodes(nidx,3);
        Alpha = PanelDat.RotAngle{i};
        AxisVec = PanelDat.HingeAxis{i}';
        RefPt = PanelDat.HingePt{i}';
        [X,Y,Z] = AxisRotating(X0,Y0,Z0,Alpha,AxisVec,RefPt);
        PanelDat.Nodes(nidx,:) = [X,Y,Z];
    end
end

% Convert Mesh format
node = PanelDat.Nodes; % 3D Node position of all aerodynamic panels (size = Number_of_nodes x 3)
ele = PanelDat.WingPanel; % Node indices of all aerodynamic panels (size = Number_of_panels x 4)

% Generate Horse Shoes Panel Data
panel_vr = mesh2panel(node,ele);
[Sc,Sm,Si,So,S,pspan,pchord,normvec] = lattice_setup(panel_vr);

%%%%%%%%% sub-functions %%%%%%%
function [iInterval,iSubSurf] = getWingSec(AC)
iInterval = [0 1];
for i = 1:size(AC.SubSurf,2)
    iInterval = [iInterval AC.SubSurf(i).SpanData(1:2)];
end
iInterval = sort(unique(iInterval));
for i = 1:(length(iInterval)-1)
    iSubSurf(i) = 0;% no flap
    for j = 1:size(AC.SubSurf,2)
        if AC.SubSurf(j).SpanData(1) == iInterval(i) && AC.SubSurf(j).SpanData(2) == iInterval(i+1)
            iSubSurf(i) = j;
            break
        end
    end
end