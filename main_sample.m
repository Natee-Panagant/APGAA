clear all; close all; clc; 

addpath(fullfile(pwd,'panel_tools'))
addpath(fullfile(pwd,'vlm_dlm'))
addpath(fullfile(pwd,'examples'))

% Name of Input function file
nameACscript = 'ex5_AC_full';

%% 
% Input function file
[AC,FC]= feval(nameACscript);

%%
% Generate Aerodynamic Panel
[PanelDat,FC,Sc,Sm,Si,So,S,pspan,pchord,normvec] = ...
                PanelGen(nameACscript); 
% Panel generation -> Input a string of input filename which is 'ex_simple_wing' in this case

% Vortex Lattice Method (VLM)
[D0,A,GAMMA,F_VLM,M_VLM] = ...
    VLM(FC.rG,FC.M,FC.Qinf,FC.rho_air,Sc,Sm,Si,So,S,pspan,normvec);

% Doublet Lattice Method (DLM)
D = DLM(Sc,Si,Sm,So,FC.M,FC.k,pchord,D0); % Calculate D matrix with DLM


%%
% Interpreting VLM calculatation
W = (FC.Qinf(1)*normvec(:,1)+...
        FC.Qinf(2)*normvec(:,2)+...
        FC.Qinf(3)*normvec(:,3))/norm(FC.Qinf); % Calculate Downwash

% Calculate Cp of steady part
Ajj_VLM = D0;
Qjj_VLM = -inv(Ajj_VLM);
Cp_VLM = Qjj_VLM*W;

% Interpreting DLM calculatation
% Calculate Cp of unsteady part
Nk = numel(FC.k); % Number of Reduced Frequencies
Cp_DLM = cell(1,Nk);

% Calculate Cp of each reduced frequency
for i = 1:Nk
    Cp_DLM{i} = -inv(D{i})*W;
end


%%
% Visualization and Monitoring
% Plotting result
plot_panel(PanelDat);
set(gca,'View',[-45,45])

plot_F(PanelDat,F_VLM);
set(gca,'View',[-45,45])

GOBJ1 = plot_Cp(PanelDat,Cp_VLM);
set(GOBJ1.ax,'View',[0,90])

GOBJ2 = plot_Cp(PanelDat,Cp_DLM,FC.k);
for i = 1:size(GOBJ2.ax,1)
    set(GOBJ2.ax{i,1},'View',[0,90])
    set(GOBJ2.ax{i,2},'View',[0,90])
end

% Display Forces and Moments summation
F_VLM_total = sum(F_VLM,1);
M_VLM_total = sum(M_VLM,1);

disp(['Fx_total = ' num2str(F_VLM_total(1),'%1.2f')]);
disp(['Fy_total = ' num2str(F_VLM_total(2),'%1.2f')]);
disp(['Fz_total = ' num2str(F_VLM_total(3),'%1.2f')]);
disp(['Mx_total = ' num2str(M_VLM_total(1),'%1.2f')]);
disp(['My_total = ' num2str(M_VLM_total(2),'%1.2f')]);
disp(['Mz_total = ' num2str(M_VLM_total(3),'%1.2f')]);

return

