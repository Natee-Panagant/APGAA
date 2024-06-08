function plot_panel(PanelDat)
clrl = get_color_list;

% Plot Aerodynamic Panels
figure(1);clf;hold on;
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
Nsurf = numel(PanelDat.SurfPanelID);
SurfName=PanelDat.SurfName;
for i=1:Nsurf
    if strcmp(SurfName{i}(end-1:end),'_R') || strcmp(SurfName{i}(end-1:end),'_L')
        SurfName{i}(end-1:end)=[];
    end
end
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

for i=1:numel(PanelDat.SurfPanelID)
    id = PanelDat.SurfPanelID{i};
    x = [node(ele(id,1),1) , node(ele(id,2),1) , node(ele(id,3),1) , node(ele(id,4),1) , node(ele(id,1),1)]';
    y = [node(ele(id,1),2) , node(ele(id,2),2) , node(ele(id,3),2) , node(ele(id,4),2) , node(ele(id,1),2)]';
    z = [node(ele(id,1),3) , node(ele(id,2),3) , node(ele(id,3),3) , node(ele(id,4),3) , node(ele(id,1),3)]';
    patch(x,y,z,zeros(size(x)),'facecolor',clrl{clr_idx(i)},'facealpha',0.5);
    if numel(PanelDat.HingePt{i})>0
        for j=1:size(PanelDat.HingePt{i},1)
            arrow3dRoundHead(PanelDat.HingePt{i}(j,:)',PanelDat.HingePt{i}(j,:)'+PanelDat.HingeAxis{i}(j,:)')
        end
    end
end
title('Aerodynamic Panels');
axis equal;
view(-30,30);
end

% Sub-Function
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