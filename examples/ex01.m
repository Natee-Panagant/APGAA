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
[PanelDat,FC,Sc,Sm,Si,So,S,pspan,pchord,normvec] = PanelGen('AcModel00'); % Panel generation -> Input a string of input filename which is 'ex_simple_wing' in this case

% Vortex Lattice Method (VLM)
[D0,A,GAMMA,F_VLM,M_VLM] = VLM(FC.rG,FC.M,FC.Qinf,FC.rho_air,Sc,Sm,Si,So,S,pspan,normvec);
W = (FC.Qinf(1)*normvec(:,1)+FC.Qinf(2)*normvec(:,2)+FC.Qinf(3)*normvec(:,3))/norm(FC.Qinf);

% Calculate Cp of steady part
Ajj_VLM=D0;
Qjj_VLM = -inv(Ajj_VLM);
Cp_VLM = Qjj_VLM*W;

% Doublet Lattice Method (DLM)
D = DLM(Sc,Si,Sm,So,FC.M,FC.k,pchord,D0); % Calculate D matrix with DLM
W = (FC.Qinf(1)*normvec(:,1)+FC.Qinf(2)*normvec(:,2)+FC.Qinf(3)*normvec(:,3))/norm(FC.Qinf); % Calculate Downwash

% Calculate Cp of unsteady part
Nk = numel(FC.k); % Number of Reduced Frequencies
Cp_DLM = cell(1,Nk);
for i = 1:Nk
    Cp_DLM{i} = -inv(D{i})*W; % Calculate Cp of each reduced frequency
end

F_VLM_total = sum(F_VLM,1);
M_VLM_total = sum(M_VLM,1);
disp(['Fx_total = ' num2str(F_VLM_total(1),'%1.2f')]);
disp(['Fy_total = ' num2str(F_VLM_total(2),'%1.2f')]);
disp(['Fz_total = ' num2str(F_VLM_total(3),'%1.2f')]);
disp(['Mx_total = ' num2str(M_VLM_total(1),'%1.2f')]);
disp(['My_total = ' num2str(M_VLM_total(2),'%1.2f')]);
disp(['Mz_total = ' num2str(M_VLM_total(3),'%1.2f')]);
% Plot results
plot_panel(PanelDat);
% plot_Cp(PanelDat,Cp_VLM);
% plot_Cp(PanelDat,Cp_DLM,State.k);