function [Sc,Sm,Si,So,S,pspan,pchord,normvec,panel0] = lattice_setup(panel0)
% Generate vortex lattice panel
npanel = size(panel0,1);
panel_vr = panel0;
Sc = zeros(npanel,3);% Sc in DLM is also called Collocation points in VLM
Sm = zeros(npanel,3);% Sm in DLM is also called Vortex points in VLM
Si = zeros(npanel,3);% Si points in DLM (A point at the left side of each bound vortex)
So = zeros(npanel,3);% So points in DLM (A point at the right side of each bound vortex)
normvec = zeros(npanel,3);
for i = 1:size(panel0,1)
    vd1 = 0.5*(panel0(i,4,:)-panel0(i,1,:));
    vd2 = 0.5*(panel0(i,3,:)-panel0(i,2,:));%->already shift from PanelGen02
    Sc(i,:) = 0.5*(panel_vr(i,1,:)+vd1+panel_vr(i,2,:)+vd2);
    Sm(i,:) = 0.5*(panel_vr(i,1,:)+panel_vr(i,2,:));
    if panel_vr(i,1,2)>panel_vr(i,2,2)
        Si(i,:) = panel_vr(i,2,:);
        So(i,:) = panel_vr(i,1,:);
    elseif panel_vr(i,1,2)<panel_vr(i,2,2)
        Si(i,:) = panel_vr(i,1,:);
        So(i,:) = panel_vr(i,2,:);
    else
        if panel_vr(i,1,3)>panel_vr(i,2,3)
            Si(i,:) = panel_vr(i,2,:);
            So(i,:) = panel_vr(i,1,:);
        elseif panel_vr(i,1,3)<panel_vr(i,2,3)
            Si(i,:) = panel_vr(i,1,:);
            So(i,:) = panel_vr(i,2,:);
        else
            error('Mesh error, Please check!');
        end
    end
    
    normvec(i,:) = cross(Sc(i,:)-So(i,:),Sc(i,:)-Si(i,:));
    normvec(i,:) = normvec(i,:)/norm(normvec(i,:));% Normal vector of each panel
end

P1 = permute(panel0(:,1,:),[1 3 2]);
P2 = permute(panel0(:,2,:),[1 3 2]);
P3 = permute(panel0(:,3,:),[1 3 2]);
P4 = permute(panel0(:,4,:),[1 3 2]);
A1 = cross(P2-P1,P4-P1);
A2 = cross(P2-P3,P4-P3);
S = 0.5*sqrt(A1(:,1).^2+A1(:,2).^2+A1(:,3).^2)+0.5*sqrt(A2(:,1).^2+A2(:,2).^2+A2(:,3).^2);% Panel area
dels = 0.5*(P1(:,2:3)+P4(:,2:3))-0.5*(P2(:,2:3)+P3(:,2:3));

pspan = sqrt(dels(:,1).^2+dels(:,2).^2);% Panel span
pchord = abs(S./pspan);% Panel chord
end