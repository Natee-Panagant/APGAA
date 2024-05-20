clear all;close all;clc;
format short g

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
[AC, PanelDat, AcControl, AcTrim]=PanelGen04('Geo_DLR_Validation_FullAC',state);

% Convert Mesh format
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
Npanel = size(ele,1);

% Generate Horse Shoes Panel Data
panel_vr = mesh2panel(node,ele);
[Sc,Sm,Si,So,S,pspan,pchord,normvec]=lattice_setup2(panel_vr);

%%%%%%%%%%%%%%%%%%
% VLM Validation %
%%%%%%%%%%%%%%%%%%
disp(repmat('%',1,60));
disp(' ');
disp('VLM Validation');
disp('inHouse-VLM_DLR vs PanelAero (DLR)');
disp(' ');

% inHouse-VLM_DLR
[D0,A,GAMMA,RHS,qxV,qyV,qzV,F_VLM]=VLM_DLR(M,Q,rho_air,Sc,Sm,Si,So,S,pspan,normvec);
wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;
Qjj_inHouse_VLM_DLR = -inv(D0);
Cp_inHouse_VLM_DLR = Qjj_inHouse_VLM_DLR*wj;


% PanelAero (DLR) - load result only
P = csvread('VLM_PanelAero_DLR_rst.csv');
Cp_PanelAero_VLM = P(:,5);
Qjj_PanelAero_VLM = csvread('VLM_PanelAero_DLR_Qjj.csv');
disp(['Total absolute diff of Cp = ' num2str(sum(sum(sum(abs(Cp_PanelAero_VLM-Cp_inHouse_VLM_DLR)))))]);
disp(['Total absolute diff of Qjj (AIC) = ' num2str(sum(sum(sum(abs(Qjj_PanelAero_VLM-Qjj_inHouse_VLM_DLR)))))]);
disp(' ');

%%%%%%%%%%%%%%%%%%
% DLM Validation %
%%%%%%%%%%%%%%%%%%
disp(repmat('%',1,60));
disp(' ');
disp('DLM Validation');
disp('inHouse-VLM_DLR vs PanelAero(DLR)');
disp(' ');

% inHouse-DLM_DLR
tic
D=Dmatrix_DLRfullvec(Sc,Si,Sm,So,M,k,normvec,pspan,pchord,D0);
%     D=Dmatrix_DLR(Sc,Si,Sm,So,M,k,normvec,D0);
toc

wj = (Q(1)*normvec(:,1)+Q(2)*normvec(:,2)+Q(3)*normvec(:,3))/Qinfabs;

Qjj_inHouse_DLM_DLR1 = -inv(D{1,1});
Cp_inHouse_DLM_DLR1 = Qjj_inHouse_DLM_DLR1*wj;

Qjj_inHouse_DLM_DLR2 = -inv(D{1,2});
Cp_inHouse_DLM_DLR2 = Qjj_inHouse_DLM_DLR2*wj;

Qjj_inHouse_DLM_DLR3 = -inv(D{1,3});
Cp_inHouse_DLM_DLR3 = Qjj_inHouse_DLM_DLR3*wj;

% PanelAero (DLR) - load result only
P = csvread('DLM_PanelAero_DLR_rst.csv');
Cp_PanelAero_DLM1 = csvread('DLM_PanelAero_DLR_Cp1_real.csv') + sqrt(-1)*csvread('DLM_PanelAero_DLR_Cp1_imag.csv');
Cp_PanelAero_DLM2 = csvread('DLM_PanelAero_DLR_Cp2_real.csv') + sqrt(-1)*csvread('DLM_PanelAero_DLR_Cp2_imag.csv');
Cp_PanelAero_DLM3 = csvread('DLM_PanelAero_DLR_Cp3_real.csv') + sqrt(-1)*csvread('DLM_PanelAero_DLR_Cp3_imag.csv');

Qjj_PanelAero_DLM1 = csvread('DLM_PanelAero_DLR_Qjj1_real.csv') + sqrt(-1)*csvread('DLM_PanelAero_DLR_Qjj1_imag.csv');
Qjj_PanelAero_DLM2 = csvread('DLM_PanelAero_DLR_Qjj2_real.csv') + sqrt(-1)*csvread('DLM_PanelAero_DLR_Qjj2_imag.csv');
Qjj_PanelAero_DLM3 = csvread('DLM_PanelAero_DLR_Qjj3_real.csv') + sqrt(-1)*csvread('DLM_PanelAero_DLR_Qjj3_imag.csv');

%
disp(' ');
disp(['Total real diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(sum(sum(sum(abs(real(Cp_PanelAero_DLM1)-real(Cp_inHouse_DLM_DLR1))))))]);
disp(['Total imag diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(sum(sum(sum(abs(imag(Cp_PanelAero_DLM1)-imag(Cp_inHouse_DLM_DLR1))))))]);
disp(['Total real diff of Cp2 (k = ' num2str(k(2)) ') = ' num2str(sum(sum(sum(abs(real(Cp_PanelAero_DLM2)-real(Cp_inHouse_DLM_DLR2))))))]);
disp(['Total imag diff of Cp2 (k = ' num2str(k(2)) ') = ' num2str(sum(sum(sum(abs(imag(Cp_PanelAero_DLM2)-imag(Cp_inHouse_DLM_DLR2))))))]);
disp(['Total real diff of Cp3 (k = ' num2str(k(3)) ') = ' num2str(sum(sum(sum(abs(real(Cp_PanelAero_DLM3)-real(Cp_inHouse_DLM_DLR3))))))]);
disp(['Total imag diff of Cp3 (k = ' num2str(k(3)) ') = ' num2str(sum(sum(sum(abs(imag(Cp_PanelAero_DLM3)-imag(Cp_inHouse_DLM_DLR3))))))]);
disp(' ');
disp(['Total real diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(sum(sum(sum(abs(real(Qjj_PanelAero_DLM1)-real(Qjj_inHouse_DLM_DLR1))))))]);
disp(['Total imag diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(sum(sum(sum(abs(imag(Qjj_PanelAero_DLM1)-imag(Qjj_inHouse_DLM_DLR1))))))]);
disp(['Total real diff of Qjj2 (k = ' num2str(k(2)) ') = ' num2str(sum(sum(sum(abs(real(Qjj_PanelAero_DLM2)-real(Qjj_inHouse_DLM_DLR2))))))]);
disp(['Total imag diff of Qjj2 (k = ' num2str(k(2)) ') = ' num2str(sum(sum(sum(abs(imag(Qjj_PanelAero_DLM2)-imag(Qjj_inHouse_DLM_DLR2))))))]);
disp(['Total real diff of Qjj3 (k = ' num2str(k(3)) ') = ' num2str(sum(sum(sum(abs(real(Qjj_PanelAero_DLM3)-real(Qjj_inHouse_DLM_DLR3))))))]);
disp(['Total imag diff of Qjj3 (k = ' num2str(k(3)) ') = ' num2str(sum(sum(sum(abs(imag(Qjj_PanelAero_DLM3)-imag(Qjj_inHouse_DLM_DLR3))))))]);
disp(' ');
disp(['%Max real diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((real(Cp_PanelAero_DLM1)-real(Cp_inHouse_DLM_DLR1)))))))]);
disp(['%Max imag diff of Cp1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((imag(Cp_PanelAero_DLM1)-imag(Cp_inHouse_DLM_DLR1)))))))]);
disp(['%Max real diff of Cp1 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((real(Cp_PanelAero_DLM2)-real(Cp_inHouse_DLM_DLR2)))))))]);
disp(['%Max imag diff of Cp1 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((imag(Cp_PanelAero_DLM2)-imag(Cp_inHouse_DLM_DLR2)))))))]);
disp(['%Max real diff of Cp1 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((real(Cp_PanelAero_DLM3)-real(Cp_inHouse_DLM_DLR3)))))))]);
disp(['%Max imag diff of Cp1 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((imag(Cp_PanelAero_DLM3)-imag(Cp_inHouse_DLM_DLR3)))))))]);
disp(' ');
disp(['%Max real diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((real(Qjj_PanelAero_DLM1)-real(Qjj_inHouse_DLM_DLR1)))))))]);
disp(['%Max imag diff of Qjj1 (k = ' num2str(k(1)) ') = ' num2str(max(max(max(abs((imag(Qjj_PanelAero_DLM1)-imag(Qjj_inHouse_DLM_DLR1)))))))]);
disp(['%Max real diff of Qjj1 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((real(Qjj_PanelAero_DLM2)-real(Qjj_inHouse_DLM_DLR2)))))))]);
disp(['%Max imag diff of Qjj1 (k = ' num2str(k(2)) ') = ' num2str(max(max(max(abs((imag(Qjj_PanelAero_DLM2)-imag(Qjj_inHouse_DLM_DLR2)))))))]);
disp(['%Max real diff of Qjj1 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((real(Qjj_PanelAero_DLM3)-real(Qjj_inHouse_DLM_DLR3)))))))]);
disp(['%Max imag diff of Qjj1 (k = ' num2str(k(3)) ') = ' num2str(max(max(max(abs((imag(Qjj_PanelAero_DLM3)-imag(Qjj_inHouse_DLM_DLR3)))))))]);
disp(' ');
disp(repmat('%',1,60));

%Visualization
atol = 1e-4;
rtol = 1e-3;

%%%%%%%%%%%%%%%%%
% dCp Diff Plot %
%%%%%%%%%%%%%%%%%
disp('Absolute (atol) + Relative (rtol) error check');
disp('If abs(Cp1-Cp2) <= atol+rtol*abs(Cp2) = converged');
disp(['atol = ' num2str(atol) ', rtol = ' num2str(rtol)]);
Idiff{1} = abs(real(Cp_PanelAero_DLM1)-real(Cp_inHouse_DLM_DLR1)) > atol+rtol*abs(real(Cp_PanelAero_DLM1));%k1 real
Idiff{2} = abs(real(Cp_PanelAero_DLM2)-real(Cp_inHouse_DLM_DLR2)) > atol+rtol*abs(real(Cp_PanelAero_DLM2));%k2 real
Idiff{3} = abs(real(Cp_PanelAero_DLM3)-real(Cp_inHouse_DLM_DLR3)) > atol+rtol*abs(real(Cp_PanelAero_DLM3));%k3 real

Idiff{4} = abs(imag(Cp_PanelAero_DLM1)-imag(Cp_inHouse_DLM_DLR1)) > atol+rtol*abs(imag(Cp_PanelAero_DLM1));%k1 imag
Idiff{5} = abs(imag(Cp_PanelAero_DLM2)-imag(Cp_inHouse_DLM_DLR2)) > atol+rtol*abs(imag(Cp_PanelAero_DLM2));%k2 imag
Idiff{6} = abs(imag(Cp_PanelAero_DLM3)-imag(Cp_inHouse_DLM_DLR3)) > atol+rtol*abs(imag(Cp_PanelAero_DLM3));%k3 imag

x = panel_vr(:,1:4,1)';
y = panel_vr(:,1:4,2)';
z = panel_vr(:,1:4,3)';

cmap = [1 1 1
    1 0 0];
colormap(cmap);

figure(1);
for i=1:6
    subplot(2,3,i);hold on;
    c=zeros(1,Npanel);
    c(Idiff{i})=1;
    Ndiff = sum(c);
    switch i
        case 1
            title(['ddCp_real (@k=0.001) (Ndiff = ' num2str(Ndiff) ')'],'interpret','none');
        case 2
            title(['ddCp_real (@k=0.6) (Ndiff = ' num2str(Ndiff) ')'],'interpret','none');
        case 3
            title(['ddCp_real (@k=1.4) (Ndiff = ' num2str(Ndiff) ')'],'interpret','none');
        case 4
            title(['ddCp_imag (@k=0.001) (Ndiff = ' num2str(Ndiff) ')'],'interpret','none');
        case 5
            title(['ddCp_imag (@k=0.6) (Ndiff = ' num2str(Ndiff) ')'],'interpret','none');
        case 6
            title(['ddCp_imag (@k=1.4) (Ndiff = ' num2str(Ndiff) ')'],'interpret','none');
    end
    
    fill3(x,y,z,c);
    set(gca,'CLim',[0 1]);
    axis equal
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-30,30);
end
set(gcf,'windowstate','maximize')
%%%%%%
return
%%%%%%