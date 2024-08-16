function [gobj,qobj] = plot_F(rst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Aerodynamic Loadings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PanelDat = rst.PanelDat;
F = rst.F_VLM;

gobji = plot_panel(rst);
P = PanelDat.VtxPt;

hold on
set(gobji.Children(contains(get(gobji.Children,'Tag'),'Arrow')),'Visible','off');
set(gobji.Children(contains(get(gobji.Children,'Tag'),'Patch')),'FaceAlpha',0);
set(get(gca,'Title'),'String','');

xyz = P; % Sc,So
uvw = F;
x = xyz(:,1); u = uvw(:,1);
y = xyz(:,2); v = uvw(:,2);
z = xyz(:,3); w = uvw(:,3);

% u = rescale(u,0.1,1);
% v = rescale(v,0.1,1);
% w = rescale(w,0.1,1);

q1 = quiver3(x,y,z,u,zeros(size(u)),zeros(size(u)),'r');
q2 = quiver3(x,y,z,zeros(size(v)),v,zeros(size(v)),'g');
q3 = quiver3(x,y,z,zeros(size(w)),zeros(size(w)),w,'b');


set(q1,'Tag','Quiver_Force_X');
set(q2,'Tag','Quiver_Force_Y');
set(q3,'Tag','Quiver_Force_Z');

qobj = [q1,q2,q3]';

lgd = legend(gca,{'Fx','Fy','Fz'});

xlabel('Xaxis (m)')
ylabel('Yaxis (m)')
zlabel('Zaxis (m)')
axis equal;
view(-45,45);
set(gcf,'windowstate','maximize');

box on
grid on

hold off

gobj = gobji;
qobj;

end