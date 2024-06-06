function plot_panel(PanelDat,varargin)
if numel(varargin)>0
    AC = varargin{1};
    cs_flag = 1;
end

clrl = get_color_list;

% Plot Aerodynamic Panels
figure(1);clf;hold on;
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
for i=1:numel(PanelDat.SurfID)
    id = PanelDat.SurfID{i};
    x = [node(ele(id,1),1) , node(ele(id,2),1) , node(ele(id,3),1) , node(ele(id,4),1) , node(ele(id,1),1)]';
    y = [node(ele(id,1),2) , node(ele(id,2),2) , node(ele(id,3),2) , node(ele(id,4),2) , node(ele(id,1),2)]';
    z = [node(ele(id,1),3) , node(ele(id,2),3) , node(ele(id,3),3) , node(ele(id,4),3) , node(ele(id,1),3)]';
    patch(x,y,z,zeros(size(x)),'facecolor',clrl{i});
end

if cs_flag ==1
    for i=1:numel(AC)
        if isfield(AC(i),'SubSurf')
            for j=1:numel(AC(i).SubSurf)
                SS=AC(i).SubSurf(j);
                arrow3dRoundHead(SS.HingePt_R',SS.HingePt_R'+SS.RotAxis_R');
                if isfield(SS,'HingePt_L')
                    arrow3dRoundHead(SS.HingePt_L',SS.HingePt_L'+SS.RotAxis_L');
                end
            end
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