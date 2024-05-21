function plot_panel(PanelDat)
figure(1);clf;hold on;
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
x = [node(ele(:,1),1) , node(ele(:,2),1) , node(ele(:,3),1) , node(ele(:,4),1) , node(ele(:,1),1)]';
y = [node(ele(:,1),2) , node(ele(:,2),2) , node(ele(:,3),2) , node(ele(:,4),2) , node(ele(:,1),2)]';
z = [node(ele(:,1),3) , node(ele(:,2),3) , node(ele(:,3),3) , node(ele(:,4),3) , node(ele(:,1),3)]';
plot3(x,y,z,'k-');
axis equal;
view(-30,30);
title('Aerodynamic Panels');
end

