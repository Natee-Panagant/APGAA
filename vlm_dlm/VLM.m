%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% developed by: Natee Panagant                       %
%      Address: Department of Mechanical Engineering %
%               Faculty of Engineering               %
%               Khon Kaen University                 %
%               Thailand                             %
%               40002                                %
%        Email: natepa@kku.ac.th                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D0,A,GAMMA,F,M] = VLM(rG,Mach,Q,rho,Sc,Sm,Si,So,S,pspan,normvec)
%%%%%%%%%%%%%%%%%%%%%%%%%
%               Vortex Lattice Method               %
%%%%%%%%%%%%%%%%%%%%%%%%%

%## Implementation
% - Voß [1] 
% - Kotikalpudi [2,3] public available VLM
%
%## Formulation and derivation of VLM
% - Katz and Plotkin [4]

beta = sqrt(1-Mach^2); % Prandtl-Glauert correction to include compressibility effect
% beta suggested by [5] Hedman

% Find Induce Velocity
A = Aassem(Sc,Si,So,normvec,beta);
[~,IV,VortexD] = Aassem(Sm,Si,So,normvec,beta);

W = -sum(Q.*normvec,2);
GAMMA = A\W;

Vx = IV(:,:,1)*GAMMA;
Vy = IV(:,:,2)*GAMMA;
Vz = IV(:,:,3)*GAMMA;
Vflow = Q+[Vx,Vy,Vz];
F = rho*cross(Vflow,VortexD,2).*repmat(GAMMA,1,3);
M = cross(Sm-rG,F,2);

%Calculate D0 for VLM convergence mode in DLM for better accuracy at zero reduced frequency (k=0)
D0 = A*diag(S./pspan/2); % For mapping Angle of Attack (AoA) to Pressure Coefficient (Cp) *-> 0.5*panel_area/panel_span (0.5*A/b) in DLR_DLM formulation
end

%%%%%%%%%%%%%%%%
% Sub-Function %
%%%%%%%%%%%%%%%%
function [A,IV,VortexD]=Aassem(C,Si,So,normvec,beta)
npanel = size(C,1);
tol = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%
%## Implementation
% - Voß [1] 
%## Formulation and derivation
% - Katz and Plotkin [4]
%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply compressibility effect (Prandtl Glauert)
Si(:,1) = Si(:,1)/beta;
So(:,1) = So(:,1)/beta;
C(:,1) = C(:,1)/beta;

% Get Leading Edge Nodes
N1 = Si;
N2 = So;

VortexD = N2-N1;

%Calculate D1 (induced velocity of bound vortex) with Biot–Savart
r0_x = repmat(N2(:,1)-N1(:,1),1,npanel)';
r0_y = repmat(N2(:,2)-N1(:,2),1,npanel)';
r0_z = repmat(N2(:,3)-N1(:,3),1,npanel)';

r1_x = C(:,1)-N1(:,1)';
r1_y = C(:,2)-N1(:,2)';
r1_z = C(:,3)-N1(:,3)';

r2_x = C(:,1)-N2(:,1)';
r2_y = C(:,2)-N2(:,2)';
r2_z = C(:,3)-N2(:,3)';

cross_r1r2_x = r1_y.*r2_z-r1_z.*r2_y;
cross_r1r2_y = r1_z.*r2_x-r1_x.*r2_z;
cross_r1r2_z = r1_x.*r2_y-r1_y.*r2_x;

abs_r1 = sqrt(r1_x.^2 + r1_y.^2 + r1_z.^2);
abs_r2 = sqrt(r2_x.^2 + r2_y.^2 + r2_z.^2);
abs_cross_r1r2_p2 = cross_r1r2_x.^2 + cross_r1r2_y.^2 + cross_r1r2_z.^2;

dot_r0r1 = r0_x.*r1_x + r0_y.*r1_y + r0_z.*r1_z;
dot_r0r2 = r0_x.*r2_x + r0_y.*r2_y + r0_z.*r2_z;

D1_x = cross_r1r2_x.*(dot_r0r1./abs_r1 - dot_r0r2./abs_r2)./(4*pi*abs_cross_r1r2_p2);
D1_y = cross_r1r2_y.*(dot_r0r1./abs_r1 - dot_r0r2./abs_r2)./(4*pi*abs_cross_r1r2_p2);
D1_z = cross_r1r2_z.*(dot_r0r1./abs_r1 - dot_r0r2./abs_r2)./(4*pi*abs_cross_r1r2_p2);

% Handle singular cases of D1
singular_idx1 = abs_r1<=tol | abs_r2<=tol | abs_cross_r1r2_p2<=tol;
D1_x(singular_idx1) = 0;
D1_y(singular_idx1) = 0;
D1_z(singular_idx1) = 0;

%Calculate D2 (induced velocity of inner horseshoe vortex)
d = sqrt(r1_y.^2+r1_z.^2);
sin_gamma = -r1_z./d;
cos_gamma = r1_y./d;
cos_beta1 = 1;
cos_beta2 = -r1_x./abs_r1;
D2_y = -sin_gamma.*(cos_beta1-cos_beta2)./(4*pi*d);
D2_z = -cos_gamma.*(cos_beta1-cos_beta2)./(4*pi*d);
% Handle singular cases of D2
singular_idx2 = abs_r1<=tol | d<=tol;
D1_y(singular_idx2) = 0;
D1_z(singular_idx2) = 0;

%Calculate D3 (induced velocity of outer horseshoe vortex)
d = sqrt(r2_y.^2+r2_z.^2);
sin_gamma = r2_z./d;
cos_gamma = -r2_y./d;
cos_beta1 = r2_x./abs_r2;
cos_beta2 = -1;
D3_y = -sin_gamma.*(cos_beta1-cos_beta2)./(4*pi*d);
D3_z = -cos_gamma.*(cos_beta1-cos_beta2)./(4*pi*d);
% Handle singular cases of D3
singular_idx3 = abs_r2<=tol | d<=tol;
D1_y(singular_idx3) = 0;
D1_z(singular_idx3) = 0;

%Calculate D (total induced velocity)
D_x = D1_x;
D_y = D1_y + D2_y + D3_y;
D_z = D1_z + D2_z + D3_z;

%Calculate A matrix
normvec_x = repmat(normvec(:,1),1,npanel);
normvec_y = repmat(normvec(:,2),1,npanel);
normvec_z = repmat(normvec(:,3),1,npanel);
A = D_x.*normvec_x + D_y.*normvec_y + D_z.*normvec_z;

%
IV = zeros(npanel,npanel,3);
IV(:,:,1) = D_x;
IV(:,:,2) = D_y;
IV(:,:,3) = D_z;
end


%%
% References 
% [1] Voß, A. (2020). An Implementation of the Vortex Lattice and the Doublet Lattice Method [Monograph]. https://github.com/DLR-AE/PanelAero
% [2] Kotikalpudi, A. (2014). Body Freedom Flutter (BFF) Doublet Lattice Method (DLM). https://hdl.handle.net/11299/165566
% [3] Kotikalpudi, A., Pfifer, H., & Balas, G. J. (2015). Unsteady Aerodynamics Modeling for a Flexible Unmanned Air Vehicle. In AIAA Atmospheric Flight Mechanics Conference (1–0). American Institute of Aeronautics and Astronautics. https://doi.org/10.2514/6.2015-2854
% [4] Katz, J., & Plotkin, A. (2001). Low-Speed Aerodynamics (2nd ed.). Cambridge University Press. https://doi.org/10.1017/CBO9780511810329
% [5] Hedman, S. G. (1966). Vortex lattice method for calculation of quasi steady state loadings on thin elastic wings in subsonic flow. Flygtekniska forsoksanstalten.




