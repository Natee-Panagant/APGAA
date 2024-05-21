clearvars;close all;clc;
%Add main path
file_dir = mfilename('fullpath');
sidx = strfind(file_dir,'\');
current_dir = file_dir(1:sidx(end));
main_dir = file_dir(1:sidx(end-1));

cd(current_dir);

addpath(genpath(main_dir));
format long g

%Generate Aerodynamic Panel
[AC, PanelDat]=PanelGen04('ex_simple_wing');

% Convert Mesh format
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
Npanel = size(ele,1);

% Generate Horse Shoes Panel Data
panel_vr = mesh2panel(node,ele);
[Sc,Sm,Si,So,S,pspan,pchord,normvec]=lattice_setup2(panel_vr);



%%%%%%%
% VLM %
%%%%%%%
[D0,A,GAMMA,RHS,qxV,qyV,qzV,F_VLM]=VLM(M,Q,rho_air,Sc,Sm,Si,So,S,pspan,normvec);
wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

Ajj_VLM=D0;
Qjj_VLM = -inv(Ajj_VLM);
Cp_VLM = Qjj_VLM*wj;

%%%%%%%
% DLM %
%%%%%%%
D=DLM(Sc,Si,Sm,So,M,k,normvec,pspan,pchord,D0);
wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

for i = 1:numel(state.k)
    Cp_DLM{i} = -inv(D{i})*wj;
end

plot_panel(PanelDat);