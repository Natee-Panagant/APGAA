function [gobj,qobj] = plot_F(rst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Aerodynamic Loadings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = rst.F_VLM;

gobji = plot_panel(rst);
P = rst.Sm;% Vortex Point (in VLM) or Sm (in DLM)

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


q1 = quiver3(x,y,z,u,zeros(size(u)),zeros(size(u)),0,'r');
q2 = quiver3(x,y,z,zeros(size(v)),v,zeros(size(v)),0,'g');
q3 = quiver3(x,y,z,zeros(size(w)),zeros(size(w)),w,0,'b');

% quiver scaling normalization
hold off
maxd = max([x;y;z]);

hU1 = get(q1,'UData');
hV1 = get(q1,'VData');
hW1 = get(q1,'WData');
hU2 = get(q2,'UData');
hV2 = get(q2,'VData');
hW2 = get(q2,'WData');
hU3 = get(q3,'UData');
hV3 = get(q3,'VData');
hW3 = get(q3,'WData');

maxf = max([hU1;hV1;hW1;hU2;hV2;hW2;hU3;hV3;hW3]);
scale = 0.5*maxd/maxf;%Scale the longest quiver as % of max axis range

set(q1,'UData',scale*hU1,'VData',scale*hV1,'WData',scale*hW1);
set(q2,'UData',scale*hU2,'VData',scale*hV2,'WData',scale*hW2);
set(q3,'UData',scale*hU3,'VData',scale*hV3,'WData',scale*hW3);
%

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