function panelDat=MappedMeshAeroNew(AC,AcSurf,varargin)
% mode = 1 -> Left side Panel index start from outside
% Left Panels || Right Panels
% 13 16 19 22    1 4 7 10
% 14 17 20 23    2 5 8 11
% 15 18 21 24    3 6 9 12

% mode = 2 -> Left side Panel index start from inside
% Left Panels || Right Panels
% 22 19 16 13    1 4 7 10
% 23 20 17 14    2 5 8 11
% 24 21 18 15    3 6 9 12

if numel(varargin)>0
    mode = varargin{1};
    if ~(mode ==1 || mode==2)
        error('incorrect panel index mode input');
    end
else
    mode = 1;
end

% Grid generation
x2=AcSurf.Right.Xp;y2=AcSurf.Right.Yp;z2=AcSurf.Right.Zp;
[nchord,nspan]=size(x2);
nspan=nspan-1;nchord=nchord-1;
nwake=3;% to be corrected
wakeLength=10;% to be modified
Symmetry=AC.Symmetry;
% panel corner points
xx1=x2(1:end-1,1:end-1);
xx2=x2(1:end-1,2:end);
xx3=x2(2:end,2:end);
xx4=x2(2:end,1:end-1);

yy1=y2(1:end-1,1:end-1);
yy2=y2(1:end-1,2:end);
yy3=y2(2:end,2:end);
yy4=y2(2:end,1:end-1);

zz1=z2(1:end-1,1:end-1);
zz2=z2(1:end-1,2:end);
zz3=z2(2:end,2:end);
zz4=z2(2:end,1:end-1);

% panel areas
npanel=nchord*nspan;
dely=reshape(abs(yy2-yy1),npanel,1);
delx=reshape(abs(0.5*(xx3+xx4)-0.5*(xx1+xx2)),npanel,1);
S=reshape(dely.*delx,npanel,1);  

% panel grid generation
nnum=cumsum(ones(nchord+1,nspan+1))+...
    repmat((nchord+1)*(0:nspan),(nchord+1),1);
nnum1=nnum(1:end-1,1:end-1);
nnum2=nnum(1:end-1,2:end);
nnum3=nnum(2:end,2:end);
nnum4=nnum(2:end,1:end-1);
WingElem=[reshape(nnum1,npanel,1) reshape(nnum2,npanel,1) ...
    reshape(nnum3,npanel,1) reshape(nnum4,npanel,1)];

% wake vortices
Xw1=x2(nchord+1,1);Yw1=y2(nchord+1,1);Zw1=z2(nchord+1,1);
Xw2=x2(nchord+1,nspan+1);Yw2=y2(nchord+1,nspan+1);Zw2=z2(nchord+1,nspan+1);
Xw3=Xw2+wakeLength;Yw3=Yw2;Zw3=Zw2;
Xw4=Xw1+wakeLength;Yw4=Yw1;Zw4=Zw1;

[sw,tw]=meshgrid(linspace(-1,1,nwake+1),linspace(-1,1,nspan+1));sw=sw';tw=tw';% wake
xw=0.25*((1-sw).*(1-tw)*Xw1+(1-sw).*(1+tw)*Xw2+(1+sw).*(1+tw)*Xw3+(1+sw).*(1-tw)*Xw4);
yw=0.25*((1-sw).*(1-tw)*Yw1+(1-sw).*(1+tw)*Yw2+(1+sw).*(1+tw)*Yw3+(1+sw).*(1-tw)*Yw4);
zw=Zw1+(Zw2-Zw1)/(Yw2-Yw1)*(yw-Yw1);

% wake panels
nnum=cumsum(ones(nwake+1,nspan+1))+...
    repmat((nwake+1)*(0:nspan),(nwake+1),1);
nnum1=nnum(1:end-1,1:end-1);
nnum2=nnum(1:end-1,2:end);
nnum3=nnum(2:end,2:end);
nnum4=nnum(2:end,1:end-1);
WakeElem=[reshape(nnum1,nwake*nspan,1) reshape(nnum2,nwake*nspan,1) ...
    reshape(nnum3,nwake*nspan,1) reshape(nnum4,nwake*nspan,1)];
WakeTimes=reshape(cumsum(ones(nchord,nspan))',nchord*nspan,1);% will be usedf for unsteady VRM
% corresponding trailing panels
TrailElem=(nchord:nchord:npanel)';

% Vortex points
nnode=(nchord+1)*(nspan+1);
xv=reshape(0.5*(xx1+xx2),npanel,1);
yv=reshape(0.5*(yy1+yy2),npanel,1);
zv=reshape(0.5*(zz1+zz2),npanel,1);

enum=cumsum(ones(nchord,nspan))+...
    repmat((nchord)*(0:nspan-1),(nchord),1);
Ring2Lift=eye(npanel);
iRing1=reshape(enum(2:end,:),(nchord-1)*nspan,1);
iRing2=reshape(enum(1:end-1,:),(nchord-1)*nspan,1);
Ring2Lift(iRing1,iRing2)=Ring2Lift(iRing1,iRing2)-eye(npanel-nspan);

% Collocation points
xc=reshape(0.25*(xx1+xx2+xx3+xx4),npanel,1);
yc=reshape(0.25*(yy1+yy2+yy3+yy4),npanel,1);
zc=reshape(0.25*(zz1+zz2+zz3+zz4),npanel,1);

% normal vectors
vec1(:,:,1)=xx3-xx1;
vec1(:,:,2)=yy3-yy1;
vec1(:,:,3)=zz3-zz1;
vec2(:,:,1)=xx2-xx4;
vec2(:,:,2)=yy2-yy4;
vec2(:,:,3)=zz2-zz4;

nv=cross(vec1,vec2,3);
magnv=sqrt(dot(nv,nv,3));
nvx=nv(:,:,1)./magnv;
nvy=nv(:,:,2)./magnv;
nvz=nv(:,:,3)./magnv;
NormV=[reshape(nvx,npanel,1) reshape(nvy,npanel,1) reshape(nvz,npanel,1)];

% output
panelDat.Nodes=[reshape(x2,nnode,1) reshape(y2,nnode,1) reshape(z2,nnode,1)];
panelDat.NodesW=[reshape(xw,(nwake+1)*(nspan+1),1) reshape(yw,(nwake+1)*(nspan+1),1) ...
    reshape(zw,(nwake+1)*(nspan+1),1)];
panelDat.VtxPt=[xv yv zv];
panelDat.Ring2Lift=Ring2Lift;
panelDat.ColPt=[xc yc zc];
panelDat.NormV=NormV;
panelDat.delx=delx;
panelDat.dely=dely;
panelDat.S=S;
panelDat.WingPanel=WingElem;
panelDat.WakePanel=WakeElem;
panelDat.TrailPanel=TrailElem;

if strcmp(Symmetry,'yes')
    panelDat.Nodes=[panelDat.Nodes
        panelDat.Nodes(:,1) -panelDat.Nodes(:,2) panelDat.Nodes(:,3)
        ];
    panelDat.NodesW=[panelDat.NodesW
        panelDat.NodesW(:,1) -panelDat.NodesW(:,2) panelDat.NodesW(:,3)
        ];
    panelDat.VtxPt=[panelDat.VtxPt
        panelDat.VtxPt(:,1) -panelDat.VtxPt(:,2) panelDat.VtxPt(:,3)
        ];
    panelDat.Ring2Lift=[Ring2Lift zeros(npanel)
                        zeros(npanel) Ring2Lift];
    panelDat.ColPt=[panelDat.ColPt
        panelDat.ColPt(:,1) -panelDat.ColPt(:,2) panelDat.ColPt(:,3)
        ];
    panelDat.NormV=[panelDat.NormV
        panelDat.NormV(:,1) -panelDat.NormV(:,2) panelDat.NormV(:,3)
        ];
    panelDat.delx=[panelDat.delx;panelDat.delx];
    panelDat.dely=[panelDat.dely;panelDat.dely];
    panelDat.S=[panelDat.S;panelDat.S];
    if mode == 1
        panelDat.WingPanel=[WingElem;nnode+WingElem];
    elseif mode ==2
        % flip the mirror panel index -> first element from outside of the left panel
        idx = reshape(fliplr(reshape(1:npanel,nchord,nspan)),[],1);
        panelDat.WingPanel=[WingElem;nnode+WingElem(idx,:)];
    end
    panelDat.WakePanel=[WakeElem;(nwake+1)*(nspan+1)+WakeElem];
    panelDat.TrailPanel=[TrailElem;npanel+TrailElem];
end
