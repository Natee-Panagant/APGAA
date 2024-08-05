function gobj = plot_panel(PanelDat,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Aerodynamic Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(varargin) > 0
    wake_flag = varargin{1};
else
    wake_flag = 0;
end
clrl = get_color_list;
figure; hold on;
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
Nsurf = numel(PanelDat.SurfPanelID);
SurfName = PanelDat.SurfName;
% Get the surface list
for i=1:Nsurf
    if strcmp(SurfName{i}(end-1:end),'_R') || strcmp(SurfName{i}(end-1:end),'_L')
        SurfName{i}(end-1:end)=[];
    end
end
% Set panel colors
c=0;
for i=1:Nsurf
    dup_idx = find(ismember(SurfName,SurfName{i})==1);
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
hg = hggroup;
for i=1:numel(PanelDat.SurfPanelID)
    %
    id = PanelDat.SurfPanelID{i};
    xyz = mesh2panel(node,ele(id,:));
    %% xxx
    % patch(xyz(:,:,1)',xyz(:,:,2)',xyz(:,:,3)',0,'facecolor',clrl{clr_idx(i)},'facealpha',0.5);
    patsurf(i,1) = patch(xyz(:,:,1)',xyz(:,:,2)',xyz(:,:,3)',0,'facecolor',clrl{clr_idx(i)},'facealpha',0.5);
    set(patsurf(i,1),'Tag',['Patch' '_' PanelDat.SurfName{i}]);
    %
    
    if numel(PanelDat.HingePt{i})>0
        for j=1:size(PanelDat.HingePt{i},1)
            %% xxx
            % arrow3dRoundHead(PanelDat.HingePt{i}(j,:)',PanelDat.HingePt{i}(j,:)'+PanelDat.HingeAxis{i}(j,:)')
            arrowobj(i,j) = arrow3dRoundHead(PanelDat.HingePt{i}(j,:)',PanelDat.HingePt{i}(j,:)'+PanelDat.HingeAxis{i}(j,:)');
            set(arrowobj(i,j),'Tag',['Arrow' '_' PanelDat.SurfName{i}]);
            %
        end
    end
    
end

if ~isempty(who('arrowobj'))

arrowobjn = [];
c=0;
    for i = 1:numel(arrowobj)
        if ~isempty(fieldnames(get(arrowobj(i))))
            c = c+1;
            arrowobjn = [arrowobjn; arrowobj(i)];
        end
    end
    arrowobj = arrowobjn;
    set(arrowobj(:),'Parent',hg);
end
set(patsurf(:),'Parent',hg);


% Plot wakes
if wake_flag == 1
    nodew = PanelDat.NodesW;
    elew = PanelDat.WakePanel;
    xyz = mesh2panel(nodew,elew);
    %% 
    % patch(xyz(:,:,1)',xyz(:,:,2)',xyz(:,:,3)',0,'facecolor',[.5 .5 .5],'facealpha',0.5);
    patwake = patch(xyz(:,:,1)',xyz(:,:,2)',xyz(:,:,3)',0,'facecolor',[.5 .5 .5],'facealpha',0.5);
    set(patwake,'Tag',['Patch' '_' 'WakePanel']);
    set(patwake(:),'Parent',hg);
end

% title('Aerodynamic Panels');
% axis equal;
% view(-30,30);
set(gcf,'windowstate','maximize');

xlabel('Xaxis (m)')
ylabel('Yaxis (m)')
zlabel('Zaxis (m)')
axis equal;
box on
grid on
hold off;

%%
gobj = hg;


%%
list_surfname = PanelDat.SurfName;
list_surfname = unique(replace(list_surfname,{'_L','_R'},{''}),'stable');

lgd_obj = [];
lgd_name = [];
list_obj = {'Patch','Arrow'};
for i = 1:numel(list_surfname)
    sel_ind = contains(get(gobj.Children,'Tag'),[list_obj{1},'_',list_surfname{i}]);
    obji = gobj.Children(sel_ind);
    if numel(obji)>1
        obji = obji(1);
    end
    lgd_obj = [lgd_obj; obji];
    lgd_name = [lgd_name; list_surfname(i)];
end
lgd = legend(lgd_obj,lgd_name);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Sub-Functions %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clrl = get_color_list()
clrl = {"#FF0000" % red
        "#00FF00" % green
        "#0000FF" % blue
        "#00FFFF" % cyan
        "#FF00FF" % magenta
        "#FFFF00" % yellow
        "#0072BD" % blue sea
        "#D95319" % orange
        "#EDB120" % yellow low light
        "#7E2F8E" % purple
        "#77AC30" % green light
        "#4DBEEE" % cyan low light
        "#A2142F"}; % red 
end