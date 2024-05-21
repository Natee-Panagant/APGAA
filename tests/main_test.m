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
clist = {'fullAC'; %1
    'flatpanel_5x5'; %2
    'flatpanel_5x5_sweep'; %3
    'flatpanel_5x5_dihedral_10'; %4
    'flatpanel_5x5_dihedral_45'; %5
    'flatpanel_5x5_vertical'; %6
    'flatpanel_5x5_dihedral_10_sweep'; %7
    'near_planar_case'; %8
    'further_away_case'}; %9
csel = 1:9;

err_idx = zeros(numel(clist),1);
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

    %%%%%%%
    % VLM %
    %%%%%%%
    disp(' ');
    disp('VLM testing - Current vs Ref_data');
    disp(' ');

    [D0,A,GAMMA,RHS,qxV,qyV,qzV,F_VLM]=VLM(M,Q,rho_air,Sc,Sm,Si,So,S,pspan,normvec);
    wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

    Ajj_VLM=D0;
    Qjj_VLM = -inv(Ajj_VLM);
    Cp_VLM = Qjj_VLM*wj;

    %%%%%%%
    % DLM %
    %%%%%%%
    disp(repmat('-',1,60));
    disp(' ');
    disp('DLM testing - Current vs Ref_data');
    disp(' ');

    D=DLM(Sc,Si,Sm,So,M,k,normvec,pspan,pchord,D0);
    wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

    Qjj_DLM1 = -inv(D{1,1});
    Cp_DLM1 = Qjj_DLM1*wj;

    Qjj_DLM2 = -inv(D{1,2});
    Cp_DLM2 = Qjj_DLM2*wj;

    Qjj_DLM3 = -inv(D{1,3});
    Cp_DLM3 = Qjj_DLM3*wj;
    
    % Load reference validated VLM+DLM results
    load(['rst_' clist{ci} '.mat']);

    % Comparison
    rst_diff(1) = max(max(max(abs(Cp_VLM_ref-Cp_VLM))));
    rst_diff(2) = max(max(max(abs(Qjj_VLM_ref-Qjj_VLM))));

    rst_diff(3) = max(max(max(abs((real(Cp_DLM1_ref)-real(Cp_DLM1))))));
    rst_diff(4) = max(max(max(abs((imag(Cp_DLM1_ref)-imag(Cp_DLM1))))));
    rst_diff(5) = max(max(max(abs((real(Cp_DLM2_ref)-real(Cp_DLM2))))));
    rst_diff(6) = max(max(max(abs((imag(Cp_DLM2_ref)-imag(Cp_DLM2))))));
    rst_diff(7) = max(max(max(abs((real(Cp_DLM3_ref)-real(Cp_DLM3))))));
    rst_diff(8) = max(max(max(abs((imag(Cp_DLM3_ref)-imag(Cp_DLM3))))));

    rst_diff(9) = max(max(max(abs((real(Qjj_DLM1_ref)-real(Qjj_DLM1))))));
    rst_diff(10) = max(max(max(abs((imag(Qjj_DLM1_ref)-imag(Qjj_DLM1))))));
    rst_diff(11) = max(max(max(abs((real(Qjj_DLM2_ref)-real(Qjj_DLM2))))));
    rst_diff(12) = max(max(max(abs((imag(Qjj_DLM2_ref)-imag(Qjj_DLM2))))));
    rst_diff(13) = max(max(max(abs((real(Qjj_DLM3_ref)-real(Qjj_DLM3))))));
    rst_diff(14) = max(max(max(abs((imag(Qjj_DLM3_ref)-imag(Qjj_DLM3))))));

    disp(['Maximum absolute diff of Cp = ' num2str(rst_diff(1))]);
    disp(['Maximum absolute diff of Qjj (AIC) = ' num2str(rst_diff(2))]);
    disp(' ');
    disp(['%Max real diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(rst_diff(3))]);
    disp(['%Max imag diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(rst_diff(4))]);
    disp(['%Max real diff of Cp1 (k = ' num2str(k(2)) ') = ' num2str(rst_diff(5))]);
    disp(['%Max imag diff of Cp1 (k = ' num2str(k(2)) ') = ' num2str(rst_diff(6))]);
    disp(['%Max real diff of Cp1 (k = ' num2str(k(3)) ') = ' num2str(rst_diff(7))]);
    disp(['%Max imag diff of Cp1 (k = ' num2str(k(3)) ') = ' num2str(rst_diff(8))]);
    disp(' ');
    disp(['%Max real diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(rst_diff(9))]);
    disp(['%Max imag diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(rst_diff(10))]);
    disp(['%Max real diff of Qjj2 (k = ' num2str(k(2)) ') = ' num2str(rst_diff(11))]);
    disp(['%Max imag diff of Qjj2 (k = ' num2str(k(2)) ') = ' num2str(rst_diff(12))]);
    disp(['%Max real diff of Qjj3 (k = ' num2str(k(3)) ') = ' num2str(rst_diff(13))]);
    disp(['%Max imag diff of Qjj3 (k = ' num2str(k(3)) ') = ' num2str(rst_diff(14))]);
    disp(' ');

    if max(rst_diff) > 1e-10
        err_idx(ci) = 1;
    end
end
disp(repmat('%',1,60));
disp(' ');
disp(repmat('*',1,60));
disp(repmat('*',1,60));
disp(repmat('*',1,60));
if sum(err_idx)>0
    disp('found error in following cases:');
    for ci=1:numel(err_idx)
        if err_idx(ci)>0
            disp(clist{ci});
        end
    end
else
    disp('All results converged, no errors found');
end
disp(repmat('*',1,60));
disp(repmat('*',1,60));
disp(repmat('*',1,60));