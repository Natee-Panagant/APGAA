function rst = run_solver(inp_filename)
%%%%%%%%%%%%%%%%
% APGAA solver %
%%%%%%%%%%%%%%%%

% Clear existed data and figures
clearvars -except inp_filename;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Generate Aerodynamic panel %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - The aerodynamic panels have to be created by PanelGen function. Inside of this function there are three main functional process including 
% - Translating structural grids into aerodynamic panels
% - Mapping and shifting the aerodynamic panels
% - Create the lattice panels from aerodynamic panels


% PanelGen
[PanelDat,Flight_Condition,Sc,Sm,Si,So,S,pspan,pchord,normvec] = PanelGen(inp_filename);

% Input variables %

% inp_filname                 (-) input filename

% Output variables %

% PanelDat          (-) Generated Panel data
% Flight_Condition  (-) Flight Condition data
% Sc                (m) [x y z] Collocation point 
% Sm                (m) [x y z] Vortex point
% Si                (m) [x y z] A point at left side of bound vortex
% So                (m) [x y z] A point at right side of bound vortex
% S                 (m^2) Area of panel
% pspan             (m)  Span of panel
% pchord            (m)  Chord of panel
% normvec           (-) [u v w] Normal vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Vortex Lattice method (VLM) analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For VLM function, the input are variables from PanelGen function where are two group of fight conditions and panel parameters analysis like vortex and collocation points, areas and spans of panel, and normal vectors.
% VLM
[D0,A,GAMMA,F_VLM,M_VLM] = VLM(Flight_Condition.rG,Flight_Condition.M,Flight_Condition.Qinf,Flight_Condition.rho_air,Sc,Sm,Si,So,S,pspan,normvec);

% Input variables %

% rG            (m) Center of mass of aircraft
% M             (-) Mach number
% Qinf          (m/s) [x y z] Velocity for static loadings
% rho_air       (kg/m^3) Density of air
% Sc            (m) [x y z] Collocation point 
% Sm            (m) [x y z] Vortex point
% Si            (m) [x y z] A point at left side of bound vortex
% So            (m) [x y z] A point at right side of bound vortex
% S             (m^2) Area of panel
% pspan         (m)  Span of panel
% normvec       (-) [u v w] Normal vector

% Output variables %

% D0            ()  Induced velocity with covering with geometrical
% A             ()  Induced velocity
% GAMMA         ()  Vortex strength
% F_VLM         (N)   [Fx Fy Fz] Force
% M_VLM         (N*m) [Mx My Mz] Moment (ref. center of mass)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Doublet Lattice method (DLM) analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - DLM function is require the input variable similar to VLM function. The reduced frequency is additive variable to solve in term of transient. The DLM function either require D0 solution from VLM function or not is necessary.
% DLM
D = DLM(Sc,Si,Sm,So,Flight_Condition.M,Flight_Condition.k,pchord,D0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Input variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M             (-) Mach number
% k             (-) Reduced frequency
% rho_air       (kg/m^3) Density of air
% Sc            (m) [x y z] Collocation point 
% Sm            (m) [x y z] Vortex point
% Si            (m) [x y z] A point at left side of bound vortex
% So            (m) [x y z] A point at right side of bound vortex
% S             (m^2) Area of panel
% pspan         (m)  Span of panel
% normvec       (-) [u v w] Normal vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Output variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D             (-) Induced velocity with covering with geometrical
%                   at each reduced frequency

%%%%%%%%%%%%%%%%%%%%%%
% 4. Result Analysis %
%%%%%%%%%%%%%%%%%%%%%%

% Calculate the downwash from normal vectors
W = (Flight_Condition.Qinf(1)*normvec(:,1) + ...
     Flight_Condition.Qinf(2)*normvec(:,2) + ...
     Flight_Condition.Qinf(3)*normvec(:,3))/norm(Flight_Condition.Qinf);

% Calculate pressure cofficients from VLM
% Calculate Cp of steady part
Ajj_VLM = D0;
Qjj_VLM = -inv(Ajj_VLM);
Cp_VLM = Qjj_VLM*W;

% Calculate pressure cofficients from DLM
% Calculate Cp of unsteady part
Nk = numel(Flight_Condition.k); % Number of Reduced Frequencies
Cp_DLM = cell(1,Nk);
for i = 1:Nk
    Cp_DLM{i} = -inv(D{i})*W; % Calculate Cp of each reduced frequency
end

% Summarizing lift (Fz), drag (Fx) and sideslip (Fy) forces and moments
% Summation of Forces and Moments
F_VLM_total = sum(F_VLM,1);
M_VLM_total = sum(M_VLM,1);

% Removed unused data
PanelDat = rmfield(PanelDat,{'NodesW','VtxPt','Ring2Lift','ColPt','NormV',...
           'delx','dely','S','WakePanel','TrailPanel','nvtx','nnode','npanel','nwake'});

% Return all results and data
rst = merge_struct;

% Display Total loads
disp(['Fx_total = ' num2str(F_VLM_total(1),'%1.2f'),' N']);
disp(['Fy_total = ' num2str(F_VLM_total(2),'%1.2f'),' N']);
disp(['Fz_total = ' num2str(F_VLM_total(3),'%1.2f'),' N']);
disp(['Mx_total = ' num2str(M_VLM_total(1),'%1.2f'),' N*m']);
disp(['My_total = ' num2str(M_VLM_total(2),'%1.2f'),' N*m']);
disp(['Mz_total = ' num2str(M_VLM_total(3),'%1.2f'),' N*m']);


% Sub-Child-Function
function rst = merge_struct
    vlist = who;
    rst = struct();
    for k = 1:numel(vlist)
        eval(['rst.' vlist{k} ' = ' vlist{k} ';']);
    end
end
end