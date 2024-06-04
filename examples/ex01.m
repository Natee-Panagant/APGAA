clearvars;close all;clc;
%Add main path
file_dir = mfilename('fullpath');
sidx = strfind(file_dir,'\');
current_dir = file_dir(1:sidx(end));
main_dir = file_dir(1:sidx(end-1));

cd(current_dir);

addpath(genpath(main_dir));
format shortEng

%Generate Aerodynamic Panel
[AC, PanelDat, State] = PanelGen04('ex_simple_wing'); % Panel generation -> Input a string of input filename which is 'ex_simple_wing' in this case
M = State.M; % Mach number
Qinf = State.Qinf; % 3D Velocity vector (size = 1 x 3)
rho_air = State.rho_air; % Air density

% Convert Mesh format
node = PanelDat.Nodes; % 3D Node position of all aerodynamic panels (size = Number_of_nodes x 3)
ele = PanelDat.WingPanel; % Node indices of all aerodynamic panels (size = Number_of_panels x 4)
Npanel = size(ele,1); % Number of panels

% Generate Horse Shoes Panel Data
panel_vr = mesh2panel(node,ele);
[Sc,Sm,Si,So,S,pspan,pchord,normvec]=lattice_setup2(panel_vr);

% Plot Aerodynamic Panels
% plot_panel(PanelDat);

% Vortex Lattice Method (VLM)
[D0,A,GAMMA,RHS,qxV,qyV,qzV,F_VLM]=VLM(State.M,State.Qinf,State.rho_air,Sc,Sm,Si,So,S,pspan,normvec);
wj = (State.Qinf(1)*normvec(:,1)+State.Qinf(2)*normvec(:,2)+State.Qinf(3)*normvec(:,3))/norm(State.Qinf);

% Calculate Cp of steady part
Ajj_VLM=D0;
Qjj_VLM = -inv(Ajj_VLM);
Cp_VLM = Qjj_VLM*wj;

% Doublet Lattice Method (DLM)
D = DLM(Sc,Si,Sm,So,State.M,State.k,normvec,pspan,pchord,D0); % Calculate D matrix with DLM
wj = (State.Qinf(1)*normvec(:,1)+State.Qinf(2)*normvec(:,2)+State.Qinf(3)*normvec(:,3))/norm(State.Qinf); % Calculate Downwash

% Calculate Cp of unsteady part
Nk = numel(State.k); % Number of Reduced Frequencies
Cp_DLM = cell(1,Nk);
for i = 1:Nk
    Cp_DLM{i} = -inv(D{i})*wj; % Calculate Cp of each reduced frequency
end

F_VLM_total = sum(F_VLM,1);
disp(['Fx_total = ' num2str(F_VLM_total(1),'%1.2f')]);
disp(['Fy_total = ' num2str(F_VLM_total(2),'%1.2f')]);
disp(['Fz_total = ' num2str(F_VLM_total(3),'%1.2f')]);

% Plot results
plot_Cp(PanelDat,Cp_VLM);
% plot_Cp(PanelDat,Cp_DLM,State.k);