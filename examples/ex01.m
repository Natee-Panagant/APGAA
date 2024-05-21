clearvars;close all;clc;
%Add main path
file_dir = mfilename('fullpath');
sidx = strfind(file_dir,'\');
current_dir = file_dir(1:sidx(end));
main_dir = file_dir(1:sidx(end-1));

cd(current_dir);

addpath(genpath(main_dir));
format long g

% Test case list


%State Parameter
state.alpha     = 5; % (deg)
state.beta      = 0; % (deg)
state.rho_air   = 1.225;  % air density (kg/m^3)
state.M         = 0.8;          % Mach number
state.CG        = [0 0 0];     % center of gravity
state.k         = [0.001 0.6 1.4]; % Nastran reduce frequencies (omega*Uinf/semichord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This inclined flow is for validation only (*In our formulation will later use -> state.Qinf = [Qinfabs 0 0])
Qinfabs = 343*state.M;
state.Qinf = Qinfabs*[cosd(state.alpha)*cosd(state.beta) -cosd(state.alpha)*sind(state.beta) sind(state.alpha)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = state.alpha;
beta = state.beta;
rho_air = state.rho_air;
M = state.M;
CG = state.CG;
k = state.k;
Q = state.Qinf;
q = 0.5*rho_air*Qinfabs^2;
%%%%%%%%%%%%%%%%%%%%%
% inHouse - VLM_DLR %
%%%%%%%%%%%%%%%%%%%%%
[AC, PanelDat]=PanelGen04('ex_simple_wing',state);

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