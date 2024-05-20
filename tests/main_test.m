clear all;close all;clc;
%Add main path
current_dir = pwd;
sidx = strfind(current_dir,'\');
main_dir = current_dir(1:sidx(end)-1);
addpath(genpath(main_dir));
format long g

% Test case list
clist = {'fullAC'; %1
    'flatpanel_5x5'; %2
    'flatpanel_5x5_sweep'; %3
    'flatpanel_5x5_dihedral_10'; %4
    'flatpanel_5x5_dihedral_45'; %5
    'flatpanel_5x5_vertical'; %6
    'flatpanel_5x5_dihedral_10_sweep'; %7
    'near_planar_case'; %8
    'further_away_case'}; %9
csel = 1;

flag_p1 = 0;%Plot Ajj
flag_p2 = 0;%Plot Qjj
flag_p3 = 0;%Plot Cp
flag_p4 = 0;%Plot VLM Loadings
for ci = csel
    cname = clist{ci};
    fpath = [pwd '\PanelAero_Validation_RST\' cname '\'];

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
    [AC, PanelDat, AcControl, AcTrim]=PanelGen04(['GEO_' cname],state);

    % Convert Mesh format
    node = PanelDat.Nodes;
    ele = PanelDat.WingPanel;
    Npanel = size(ele,1);

    % Generate Horse Shoes Panel Data
    panel_vr = mesh2panel(node,ele);
    [Sc,Sm,Si,So,S,pspan,pchord,normvec]=lattice_setup2(panel_vr);

    disp(repmat('%',1,60));
    disp(' ');
    disp(['CASE - ' cname ' Validation']);
    disp(' ');
    disp(repmat('-',1,60));

    %%%%%%%%%%%%%%%%%%
    % VLM Validation %
    %%%%%%%%%%%%%%%%%%
    disp(' ');
    disp('VLM testing - Current vs Ref_data');
    disp(' ');

    % Current version of VLM code
    [D0,A,GAMMA,RHS,qxV,qyV,qzV,F_VLM]=VLM_DLR(M,Q,rho_air,Sc,Sm,Si,So,S,pspan,normvec);
    wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

    Ajj_VLM=D0;
    Qjj_VLM = -inv(Ajj_VLM);
    Cp_VLM = Qjj_VLM*wj;

    % Load reference validated results
    load(['rst_' clist{ci} '.mat']);

    % Comparison
    disp(['Maximum absolute diff of Cp = ' num2str(max(max(max(abs(Cp_VLM_ref-Cp_VLM)))))]);
    disp(['Maximum absolute diff of Qjj (AIC) = ' num2str(max(max(max(abs(Qjj_VLM_ref-Qjj_VLM)))))]);
    disp(' ');

    %%%%%%%%%%%%%%%%%%
    % DLM Validation %
    %%%%%%%%%%%%%%%%%%
    disp(repmat('-',1,60));
    disp(' ');
    disp('DLM testing - Current vs Ref_data');
    disp(' ');

    % Current code execution
    D=Dmatrix_DLRfullvec(Sc,Si,Sm,So,M,k,normvec,pspan,pchord,D0);
    wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

    Qjj_DLM1 = -inv(D{1,1});
    Cp_DLM1 = Qjj_DLM1*wj;

    Qjj_DLM2 = -inv(D{1,2});
    Cp_DLM2 = Qjj_DLM2*wj;

    Qjj_DLM3 = -inv(D{1,3});
    Cp_DLM3 = Qjj_DLM3*wj;

    % PanelAero (DLR) - load result only
    disp(['%Max real diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((real(Cp_DLM1_ref)-real(Cp_DLM1)))))))]);
    disp(['%Max imag diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((imag(Cp_DLM1_ref)-imag(Cp_DLM1)))))))]);
    disp(['%Max real diff of Cp1 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((real(Cp_DLM2_ref)-real(Cp_DLM2)))))))]);
    disp(['%Max imag diff of Cp1 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((imag(Cp_DLM2_ref)-imag(Cp_DLM2)))))))]);
    disp(['%Max real diff of Cp1 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((real(Cp_DLM3_ref)-real(Cp_DLM3)))))))]);
    disp(['%Max imag diff of Cp1 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((imag(Cp_DLM3_ref)-imag(Cp_DLM3)))))))]);
    disp(' ');
    disp(['%Max real diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((real(Qjj_DLM1_ref)-real(Qjj_DLM1)))))))]);
    disp(['%Max imag diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((imag(Qjj_DLM1_ref)-imag(Qjj_DLM1)))))))]);
    disp(['%Max real diff of Qjj2 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((real(Qjj_DLM2_ref)-real(Qjj_DLM2)))))))]);
    disp(['%Max imag diff of Qjj2 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((imag(Qjj_DLM2_ref)-imag(Qjj_DLM2)))))))]);
    disp(['%Max real diff of Qjj3 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((real(Qjj_DLM3_ref)-real(Qjj_DLM3)))))))]);
    disp(['%Max imag diff of Qjj3 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((imag(Qjj_DLM3_ref)-imag(Qjj_DLM3)))))))]);
    disp(' ');
end
disp(repmat('%',1,60));
    