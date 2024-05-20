function Dij=Dmatrix_DLRfullvec(Pc,Pi,Pm,Po,M,k,N,pspan,pchord,varargin)
%Using I1/I2 Integral formulation from
%DLR-IB-AE-GO-2020-137
%An Implementation of the Vortex Lattice and the Doublt Lattice Method Version 1.04
%Arne VoB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dmatrix calculation with Doublet Lattice Method (DLM) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code support both planar and non-planar panel.
%
%
%       Pc = 3D collocation points, the 3/4 chord as displayed below
%            (Array -> Size = Number_of_Panels * 3)
% Pi,Pm,Po = 3D doublet points, the 1/4 of chord as displayed below
%            (Array -> Size = Number_of_Panels * 3)
%        M = Mach number (Scalar)
%        b = Semi chord length (Scalar)
%       k0 = Nastran's reduced frequency (k=omega*b/Uinf)
%            ((Vector Size = 1 * Number_of_reduced_frequency)
%        k = Reduced frequency (k=omega/Uinf)
%            ((Vector Size = 1 * Number_of_reduced_frequency)
%       D0 = D0 matrix from Vortex Lattice Method (optional)
%            (Array -> Size = Number_of_Panels * Number_of_Panels)
%
%               __   __   __   __
%              |  Po             |
%              |                 |
% Uinf --->    |  Pm        Pc   |
%              |                 |
%  ^y-axis     |__Pi __   __   __|
%  |           0  1/4  1/2  3/4  1
%  |-->x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% developed by: Natee Panagant                       %
%      Address: Department of Mechanical Engineering %
%               Faculty of Engineering               %
%               Khon Kaen University                 %
%               Thailand                             %
%               40002                                %
%        Email: natepa@kku.ac.th                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(varargin{1})>0
    D0 = varargin{1};%Include VLM convergence
    VLM_conv_flag = true;
else
    VLM_conv_flag = false;
end
% Preparing array data for vectorized computation
Np = size(Pc,1);    % Number of panels
Nk = numel(k);     % Number of reduced frequencies
NM = numel(M);      % Number of Mach numbers
beta2 = 1-M^2;
delx = repmat(pchord',Np,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *Receiving panel data is equal in each row  %
% *Sending panel data is equal in each column %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = repmat(0.5*sqrt(((Po(:,3)-Pi(:,3)).^2+(Po(:,2)-Pi(:,2)).^2))',Np,1);
e2 = e.^2;
e3 = e.^3;
e4 = e.^4;
e_2 = 0.5*e;
e0 = zeros(size(e));

xsr = Pc(:,1) - Pm(:,1)';
ysr = Pc(:,2) - Pm(:,2)';
zsr = Pc(:,3) - Pm(:,3)';

sin_gamma = (Po(:,3)-Pi(:,3))'./(2*e); %gamma is dihedral angle of each panel
cos_gamma = (Po(:,2)-Pi(:,2))'./(2*e);
tan_lambda = (Po(:,1)-Pi(:,1))'./(2*e);%lambda is sweep angle of each panel (%NeoCASS will fix the sweep angle variation dues to dihedral with different equation)
gamma = asin(sin_gamma);
gamma_sr = gamma'-gamma;

xbar = xsr;
ybar = ysr.*cos_gamma + zsr.*sin_gamma;
zbar = zsr.*cos_gamma - ysr.*sin_gamma;

ybar2 = ybar.^2;
ybar4 = ybar.^4;
zbar2 = zbar.^2;
zbar4 = zbar.^4;

tol0 = (ybar2+zbar2-e2);
ratio=2*e.*abs(zbar)./tol0;
Ip=abs(zbar) <= e*0.001;    % Planar panel index
In=~Ip;                     % Non-Planar panel index
Ia=abs(ratio) <= 0.3 & In;  % Co-Planar panel index
Ir=abs(ratio) > 0.3 & In;   % The rest (further away) panel index

Dij = cell(NM,Nk);
for i = 1:NM
    for j = 1:Nk
        if k(j)==0 && VLM_conv_flag
            Dij{i,j} = D0;
        else
            del1 = zeros(Np);
            del2 = zeros(Np);
            IND1 = tol0 > 0;
            del1(IND1) = 1;
            del2(IND1) = 0;
            IND2 = tol0 == 0;
            del1(IND2) = 0;
            del2(IND2) = 0.5;
            IND3=tol0 < 0;
            del1(IND3) = 1;
            del2(IND3) = 1;
            
            %Calculate Epsilon
            epsilon=zeros(Np);
            epsilon(Ip) = 2*e(Ip)./(ybar2(Ip)-e2(Ip));% For planar panels
            ind = 2:7;
            epsilon(Ia) = 4*e4(Ia)./tol0(Ia).^2.*...
                ((-1)^ind(1)/(2*ind(1)-1)*ratio(Ia).^(2*ind(1)-4)+...
                (-1)^ind(2)/(2*ind(2)-1)*ratio(Ia).^(2*ind(2)-4)+...
                (-1)^ind(3)/(2*ind(3)-1)*ratio(Ia).^(2*ind(3)-4)+...
                (-1)^ind(4)/(2*ind(4)-1)*ratio(Ia).^(2*ind(4)-4)+...
                (-1)^ind(5)/(2*ind(5)-1)*ratio(Ia).^(2*ind(5)-4)+...
                (-1)^ind(6)/(2*ind(6)-1)*ratio(Ia).^(2*ind(6)-4));% For co-planar panels
            ATAN = atan(ratio(Ir));
            ATAN(ATAN<0) = ATAN(ATAN<0)+pi;
            epsilon(Ir) = (e2(Ir)./zbar2(Ir)).*(1-1./ratio(Ir)).*ATAN;% For further-away panels
            
            % Calculate F
            F = zeros(Np);
            % For planar panels
            F(Ip) = del1(Ip).*epsilon(Ip);
            
            % For non-planar panels (In = Ia & Ir)
            F(In) = del1(In)*2.*e(In)./(ybar2(In)+zbar2(In)-e2(In)).*(1-epsilon(In).*(zbar2(In))./(e2(In)))+del2(In)*pi./abs(zbar(In));
            
            % Calculate Kernel function with Desmarais approximation
            [P1o,P2o]     = Q1Q2kernel( e   ,xbar,ybar,zbar,zbar2,beta2,gamma_sr,tan_lambda,M(i),k(j),VLM_conv_flag);
            [P1oh,P2oh]   = Q1Q2kernel( e_2 ,xbar,ybar,zbar,zbar2,beta2,gamma_sr,tan_lambda,M(i),k(j),VLM_conv_flag);
            [P10,P20]     = Q1Q2kernel( e0   ,xbar,ybar,zbar,zbar2,beta2,gamma_sr,tan_lambda,M(i),k(j),VLM_conv_flag);
            [P1ih,P2ih]   = Q1Q2kernel(-e_2 ,xbar,ybar,zbar,zbar2,beta2,gamma_sr,tan_lambda,M(i),k(j),VLM_conv_flag);
            [P1i,P2i]     = Q1Q2kernel(-e   ,xbar,ybar,zbar,zbar2,beta2,gamma_sr,tan_lambda,M(i),k(j),VLM_conv_flag);
            
            
            
            % Calculate Quartic Coefficients (A1,B1,C1,D1,E1) for Q1
            A1=(8*P1ih)./(3*e2) - P1i./(6*e2) - (5*P10)./e2 - P1o./(6*e2) + (8*P1oh)./(3*e2);
            B1=P1i./(6.*e) - (4*P1ih)./(3*e) - P1o./(6*e) + (4*P1oh)./(3*e);
            C1=P10;
            D1=(4*P1ih)./(3*e3) - (2*P1i)./(3*e3) + (2*P1o)./(3*e3) - (4*P1oh)./(3*e3);
            E1=(4*P10)./e4 + (2*P1i)./(3*e4) - (8*P1ih)./(3*e4) + (2*P1o)./(3*e4) - (8*P1oh)./(3*e4);
            
            % Calculate Quartic Coefficients (A2,B2,C2,D2,E2) for Q2
            A2=(8*P2ih)./(3*e2) - P2i./(6*e2) - (5*P20)./e2 - P2o./(6*e2) + (8*P2oh)./(3*e2);
            B2=P2i./(6*e) - (4*P2ih)./(3*e) - P2o./(6*e) + (4*P2oh)./(3*e);
            C2=P20;
            D2=(4*P2ih)./(3*e3) - (2*P2i)./(3*e3) + (2*P2o)./(3*e3) - (4*P2oh)./(3*e3);
            E2=(4*P20)./e4 + (2*P2i)./(3*e4) - (8*P2ih)./(3*e4) + (2*P2o)./(3*e4) - (8*P2oh)./(3*e4);
            
            % Calculate D1ij
            termD1=(ybar2-zbar2).*A1+ybar.*B1+C1+ybar.*(ybar2-3*zbar2).*D1+(ybar4-6*ybar2.*zbar2+zbar4).*E1;
            termD2=ybar.*A1+0.5*B1+0.5*(3*ybar2-zbar2).*D1+2*ybar.*(ybar2-zbar2).*E1;
            termD3=A1+2*ybar.*D1+(3*ybar2-zbar2+e2/3).*E1;
            
            D1ij=delx/8/pi.*(termD1.*F+termD2.*log(((ybar-e).^2+zbar2)./((ybar+e).^2+zbar2))+2*e.*termD3);
            
            % Calculate D2_Ib
            D2ij=zeros(Np);
            Ib=1./ratio<=0.1 & In;
            term21=(delx(Ib)/16/pi)./zbar2(Ib);
            term22=((ybar2(Ib)+zbar2(Ib)).*A2(Ib)+ybar(Ib).*B2(Ib)+C2(Ib)+ybar(Ib).*(ybar2(Ib)+3*zbar2(Ib)).*D2(Ib)+(ybar4(Ib)+6*ybar2(Ib).*zbar2(Ib)-3*zbar4(Ib)).*E2(Ib)).*F(Ib);
            term23=1./((ybar(Ib)+e(Ib)).^2+zbar2(Ib));
            term24=((ybar2(Ib)+zbar2(Ib)).*ybar(Ib)+(ybar2(Ib)-zbar2(Ib)).*e(Ib)).*A2(Ib)+(ybar2(Ib)+zbar2(Ib)+ybar(Ib).*e(Ib)).*B2(Ib)+(ybar(Ib)+e(Ib)).*C2(Ib);
            term25=(ybar4(Ib)-zbar4(Ib)+(ybar2(Ib)-3*zbar2(Ib)).*ybar(Ib).*e(Ib)).*D2(Ib)+((ybar4(Ib)-2*ybar2(Ib).*zbar2(Ib)-3*zbar4(Ib)).*ybar(Ib)+ ...
                (ybar4(Ib)-6*ybar2(Ib).*zbar2(Ib)+zbar4(Ib)).*e(Ib)).*E2(Ib);
            term26=1./((ybar(Ib)-e(Ib)).^2+zbar2(Ib));
            term27=((ybar2(Ib)+zbar2(Ib)).*ybar(Ib)-(ybar2(Ib)-zbar2(Ib)).*e(Ib)).*A2(Ib)+(ybar2(Ib)+zbar2(Ib)-ybar(Ib).*e(Ib)).*B2(Ib)+(ybar(Ib)-e(Ib)).*C2(Ib);
            term28=(ybar4(Ib)-zbar4(Ib)-(ybar2(Ib)-3*zbar2(Ib)).*ybar(Ib).*e(Ib)).*D2(Ib)+((ybar4(Ib)-2*ybar2(Ib).*zbar2(Ib)-3*zbar4(Ib)).*ybar(Ib)-...
                (ybar4(Ib)-6*ybar2(Ib).*zbar2(Ib)+zbar4(Ib)).*e(Ib)).*E2(Ib);
            term29=zbar2(Ib).*log(((ybar(Ib)-e(Ib)).^2+zbar2(Ib))./((ybar(Ib)+e(Ib)).^2+zbar2(Ib))).*D2(Ib)+...
                4*zbar2(Ib).*(e(Ib)+ybar(Ib).*log(((ybar(Ib)-e(Ib)).^2+zbar2(Ib))./((ybar(Ib)+e(Ib)).^2+zbar2(Ib)))).*E2(Ib);
            
            D2ij(Ib)=term21.*(term22.*F(Ib)+term23.*(term24+term25)-term26.*(term27+term28)+term29);
            
            % Calculate D2_Ic
            Ic=1./ratio>0.1 & In;
            delta=(e(Ic)./abs(zbar(Ic))).^2.*(1-del1(Ic)-del2(Ic)*pi./ratio(Ic));
            term21=e(Ic).*delx(Ic)./(8*pi*tol0(Ic));
            term22=1./((ybar(Ic)+e(Ic)).^2+zbar2(Ic))./((ybar(Ic)-e(Ic)).^2+zbar2(Ic));
            term23=2*(ybar2(Ic)+zbar2(Ic)+e2(Ic)).*(e2(Ic).*A2(Ic)+C2(Ic));
            term24=4*ybar(Ic).*e2(Ic).*B2(Ic)+2*ybar(Ic).*(ybar4(Ic)-2*e2(Ic).*ybar2(Ic)+2*ybar2(Ic).*zbar2(Ic)+3*e4(Ic)+2*e2(Ic).*zbar2(Ic)+zbar4(Ic)).*D2(Ic);
            term25=2*(3*ybar(Ic).^6-7*e2(Ic).*ybar4(Ic)+5*ybar4(Ic).*zbar2(Ic)+6*e4(Ic).*ybar2(Ic)+6*e2(Ic).*ybar2(Ic).*zbar2(Ic) ...
                -3*e2(Ic).*zbar4(Ic)-zbar(Ic).^6+ybar2(Ic).*zbar4(Ic)-2*e4(Ic).*zbar2(Ic)).*E2(Ic);
            term26=(del1(Ic).*epsilon(Ic)+delta)./e2(Ic);
            term27=(ybar2(Ic)+zbar2(Ic)).*A2(Ic)+ybar(Ic).*B2(Ic)+C2(Ic)+ybar(Ic).*(ybar2(Ic)+3*zbar2(Ic)).*D2(Ic)+(ybar4(Ic)+6*ybar2(Ic).*zbar2(Ic)-3*zbar4(Ic)).*E2(Ic);
            term28=delx(Ic)./(8*pi).*(D2(Ic)./2.*log(((ybar(Ic)-e(Ic)).^2+zbar2(Ic))./((ybar(Ic)+e(Ic)).^2+zbar2(Ic)))+...
                2*(e(Ic)+ybar(Ic).*log(((ybar(Ic)-e(Ic)).^2+zbar2(Ic))./((ybar(Ic)+e(Ic)).^2+zbar2(Ic)))).*E2(Ic));
            
            D2ij(Ic)=term21.*(term22.*(term23+term24+term25)-term26.*term27)+term28;
            
            if VLM_conv_flag % Add VLM convergence
                Dij{i,j} = D0 + D1ij + D2ij;
            else % Ignore VLM convergence
                Dij{i,j} = D1ij + D2ij;
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%
% Sub Function %
%%%%%%%%%%%%%%%%
function [P1,P2] = Q1Q2kernel(ebar,xbar,ybar,zbar,zbar2,beta2,gamma_sr,tan_lambda,M,k,VLM_conv_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel function Computaion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Np]=size(xbar,1);
I1=zeros(Np);
I2=zeros(Np);

r1=sqrt(((ybar-ebar).^2+zbar2));
k1=k.*r1;
R=sqrt(((xbar-ebar.*tan_lambda).^2+beta2*r1.^2));
u1=(M*R-xbar+ebar.*tan_lambda)./r1/beta2;

% CASE I
IND1 = u1 >= 0;
[I1(IND1),I2(IND1)] = I_int_Desmarais(u1(IND1),k1(IND1));

% CASE II
IND2 = u1 < 0;
u0 = zeros(size(u1(IND2)));
[I10,I20] = I_int_Desmarais(u0,k1(IND2));
[I1n,I2n] = I_int_Desmarais(-u1(IND2),k1(IND2));
I1(IND2)=2*real(I10)-real(I1n)+1i*imag(I1n);
I2(IND2)=2*real(I20)-real(I2n)+1i*imag(I2n);

% Calculate K1 K2 T1 T2
K1=-I1-(M*r1./R).*exp(-1i*k1.*u1)./sqrt(1+u1.^2);
K2=3*I2+1i*k1.*M^2.*r1.^2./R.^2.*exp(-1i*k1.*u1)./sqrt(1+u1.^2)+ ...
    M*r1./R.*((1+u1.^2).*(beta2*r1.^2./R.^2)+2+M*r1.*u1./R).*...
    exp(-1i*k1.*u1)./(1+u1.^2).^1.5;

T1=cos(gamma_sr);
T2star=zbar.*(zbar.*cos(gamma_sr)+(ybar-ebar).*sin(gamma_sr));

% Handle singularity cases
IND_r0xp = r1 == 0 & xbar >= 0;
IND_r0xn = r1 == 0 & xbar < 0;
K1(IND_r0xp) = -2;
K2(IND_r0xp) = 4;
K1(IND_r0xn) = 0;
K2(IND_r0xn) = 0;

% Calculate K1 Q1 K2 Q2
if VLM_conv_flag % VLM convergence will be added later
    K10=-1-(xbar-ebar.*tan_lambda)./R;
    K20=2+((xbar-ebar.*tan_lambda)./R).*(2+beta2*r1.^2./R.^2);
    
    P1=-(K1.*exp(-1i*k.*(xbar-ebar.*tan_lambda))-K10).*T1;
    P2=-(K2.*exp(-1i*k.*(xbar-ebar.*tan_lambda))-K20).*T2star;
else % Ignore VLM convergence
    P1=-(K1.*exp(-1i*k.*(xbar-ebar.*tan_lambda))).*T1;
    P2=-(K2.*exp(-1i*k.*(xbar-ebar.*tan_lambda))).*T2star;
end
end


function [I1,I2] = I_int_Desmarais(u1,k1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desmarais Approximation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=0.009054814793;
a=[0.000319759140
    -0.000055461471
    0.002726074362
    0.005749551566
    0.031455895072
    0.106031126212
    0.406838011567
    0.798112357155
    -0.417749229098
    0.077480713894
    -0.012677284771
    0.001787032960];

im=(-1)^0.5;
expimk1u1=exp(-im*k1.*u1);
n=(1:12)';
m=1;
n_m=n/m;

I0 =  a(1).*exp(-2.^n_m(1)*b*u1).*(2.^n_m(1)*b-im*k1)./((2.^n_m(1)).^2*b^2+k1.^2)+...
    a(2).*exp(-2.^n_m(2)*b*u1).*(2.^n_m(2)*b-im*k1)./((2.^n_m(2)).^2*b^2+k1.^2)+...
    a(3).*exp(-2.^n_m(3)*b*u1).*(2.^n_m(3)*b-im*k1)./((2.^n_m(3)).^2*b^2+k1.^2)+...
    a(4).*exp(-2.^n_m(4)*b*u1).*(2.^n_m(4)*b-im*k1)./((2.^n_m(4)).^2*b^2+k1.^2)+...
    a(5).*exp(-2.^n_m(5)*b*u1).*(2.^n_m(5)*b-im*k1)./((2.^n_m(5)).^2*b^2+k1.^2)+...
    a(6).*exp(-2.^n_m(6)*b*u1).*(2.^n_m(6)*b-im*k1)./((2.^n_m(6)).^2*b^2+k1.^2)+...
    a(7).*exp(-2.^n_m(7)*b*u1).*(2.^n_m(7)*b-im*k1)./((2.^n_m(7)).^2*b^2+k1.^2)+...
    a(8).*exp(-2.^n_m(8)*b*u1).*(2.^n_m(8)*b-im*k1)./((2.^n_m(8)).^2*b^2+k1.^2)+...
    a(9).*exp(-2.^n_m(9)*b*u1).*(2.^n_m(9)*b-im*k1)./((2.^n_m(9)).^2*b^2+k1.^2)+...
    a(10).*exp(-2.^n_m(10)*b*u1).*(2.^n_m(10)*b-im*k1)./((2.^n_m(10)).^2*b^2+k1.^2)+...
    a(11).*exp(-2.^n_m(11)*b*u1).*(2.^n_m(11)*b-im*k1)./((2.^n_m(11)).^2*b^2+k1.^2)+...
    a(12).*exp(-2.^n_m(12)*b*u1).*(2.^n_m(12)*b-im*k1)./((2.^n_m(12)).^2*b^2+k1.^2);

J0 = a(1).*exp(-2.^n_m(1)*b*u1)./((2.^n_m(1)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(1)).^2*b^2-k1.^2+2.^n_m(1)*b.*u1.*((2.^n_m(1)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(1)*b+u1.*((2.^n_m(1)).^2*b^2+k1.^2)))+...%end of 1st term
    a(2).*exp(-2.^n_m(2)*b*u1)./((2.^n_m(2)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(2)).^2*b^2-k1.^2+2.^n_m(2)*b.*u1.*((2.^n_m(2)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(2)*b+u1.*((2.^n_m(2)).^2*b^2+k1.^2)))+...%end of 2nd term
    a(3).*exp(-2.^n_m(3)*b*u1)./((2.^n_m(3)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(3)).^2*b^2-k1.^2+2.^n_m(3)*b.*u1.*((2.^n_m(3)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(3)*b+u1.*((2.^n_m(3)).^2*b^2+k1.^2)))+...%end of 3rd term
    a(4).*exp(-2.^n_m(4)*b*u1)./((2.^n_m(4)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(4)).^2*b^2-k1.^2+2.^n_m(4)*b.*u1.*((2.^n_m(4)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(4)*b+u1.*((2.^n_m(4)).^2*b^2+k1.^2)))+...%end of 4th term
    a(5).*exp(-2.^n_m(5)*b*u1)./((2.^n_m(5)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(5)).^2*b^2-k1.^2+2.^n_m(5)*b.*u1.*((2.^n_m(5)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(5)*b+u1.*((2.^n_m(5)).^2*b^2+k1.^2)))+...%end of 5th term
    a(6).*exp(-2.^n_m(6)*b*u1)./((2.^n_m(6)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(6)).^2*b^2-k1.^2+2.^n_m(6)*b.*u1.*((2.^n_m(6)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(6)*b+u1.*((2.^n_m(6)).^2*b^2+k1.^2)))+...%end of 6th term
    a(7).*exp(-2.^n_m(7)*b*u1)./((2.^n_m(7)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(7)).^2*b^2-k1.^2+2.^n_m(7)*b.*u1.*((2.^n_m(7)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(7)*b+u1.*((2.^n_m(7)).^2*b^2+k1.^2)))+...%end of 7th term
    a(8).*exp(-2.^n_m(8)*b*u1)./((2.^n_m(8)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(8)).^2*b^2-k1.^2+2.^n_m(8)*b.*u1.*((2.^n_m(8)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(8)*b+u1.*((2.^n_m(8)).^2*b^2+k1.^2)))+...%end of 8th term
    a(9).*exp(-2.^n_m(9)*b*u1)./((2.^n_m(9)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(9)).^2*b^2-k1.^2+2.^n_m(9)*b.*u1.*((2.^n_m(9)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(9)*b+u1.*((2.^n_m(9)).^2*b^2+k1.^2)))+...%end of 9th term
    a(10).*exp(-2.^n_m(10)*b*u1)./((2.^n_m(10)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(10)).^2*b^2-k1.^2+2.^n_m(10)*b.*u1.*((2.^n_m(10)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(10)*b+u1.*((2.^n_m(10)).^2*b^2+k1.^2)))+...%end of 10th term
    a(11).*exp(-2.^n_m(11)*b*u1)./((2.^n_m(11)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(11)).^2*b^2-k1.^2+2.^n_m(11)*b.*u1.*((2.^n_m(11)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(11)*b+u1.*((2.^n_m(11)).^2*b^2+k1.^2)))+...%end of 11th term
    a(12).*exp(-2.^n_m(12)*b*u1)./((2.^n_m(12)).^2*b^2+k1.^2).^2.*...
    ((2.^n_m(12)).^2*b^2-k1.^2+2.^n_m(12)*b.*u1.*((2.^n_m(12)).^2*b^2+k1.^2) ...
    -im*k1.*(2*2.^n_m(12)*b+u1.*((2.^n_m(12)).^2*b^2+k1.^2)));%end of 12th term

I1 = (1-u1./(1+u1.^2).^0.5-im*k1.*I0).*expimk1u1;
I2 = 1/3*expimk1u1.*((2+im*k1.*u1).*(1-u1./(1+u1.^2).^0.5)- ...
    u1./(1+u1.^2).^1.5-im*k1.*I0+k1.^2.*J0);
end