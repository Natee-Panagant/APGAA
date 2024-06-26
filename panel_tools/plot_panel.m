function plot_panel(PanelDat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Aerodynamic Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(varargin) > 0
    wake_flag = varargin{1};
else
    wake_flag = 0;
end
clrl = get_color_list;
figure(1);clf;hold on;
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
Nsurf = numel(PanelDat.SurfPanelID);
SurfName=PanelDat.SurfName;
% Get the surface list
for i=1:Nsurf
    if strcmp(SurfName{i}(end-1:end),'_R') || strcmp(SurfName{i}(end-1:end),'_L')
        SurfName{i}(end-1:end)=[];
    end
end
% Set panel colors
c=0;
for i=1:Nsurf
    dup_idx = find( ismember(SurfName,SurfName{i})==1);
    if numel(dup_idx)>0
        if dup_idx(1)<i
            clr_idx(i)=clr_idx(dup_idx(1));
        else
            c=c+1;
            clr_idx(i)=c;
        end
    else
        c=c+1;
        clr_idx(i)=c;
    end
end
% Plot each part with different color
for i=1:numel(PanelDat.SurfPanelID)
    id = PanelDat.SurfPanelID{i};
    xyz = mesh2panel(node,ele(id,:));
    patch(xyz(:,:,1)',xyz(:,:,2)',xyz(:,:,3)',0,'facecolor',clrl{clr_idx(i)},'facealpha',0.5);
    if numel(PanelDat.HingePt{i})>0
        for j=1:size(PanelDat.HingePt{i},1)
            arrow3dRoundHead(PanelDat.HingePt{i}(j,:)',PanelDat.HingePt{i}(j,:)'+PanelDat.HingeAxis{i}(j,:)')
        end
    end
end
% Plot wakes
if wake_flag == 1
    nodew = PanelDat.NodesW;
    elew = PanelDat.WakePanel;
    xyz = mesh2panel(nodew,elew);
    patch(xyz(:,:,1)',xyz(:,:,2)',xyz(:,:,3)',0,'facecolor',[.5 .5 .5],'facealpha',0.5);
end
title('Aerodynamic Panels');
axis equal;
view(-30,30);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Sub-Functions %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clrl = get_color_list()
clrl = {"#FF0000"
        "#00FF00"
        "#0000FF"
        "#00FFFF"
        "#FF00FF"
        "#FFFF00"
        "#0072BD"
        "#D95319"
        "#EDB120"
        "#7E2F8E"
        "#77AC30"
        "#4DBEEE"
        "#A2142F"};
end