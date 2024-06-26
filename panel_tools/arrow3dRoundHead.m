function arrow3dRoundHead(p1,p2,headColor,lineColor)
% plot a 3D arrow pointing from p1 to p2
%p1=p;p2=p+dp;
if nargin==2
    headColor='r';
    lineColor='k';
elseif nargin==3
    lineColor='k';
end

x21=p1-p2;
x22=p2-[p2(1);p2(2)+1;p2(3)];
xx=cross(x21,x22);mxx=sqrt(sum(xx.*xx));
if mxx < 1e-4
    x22=p2-[p2(1);p2(2);p2(3)+1];
end
mx21=sqrt(sum(x21.*x21));
ex=x21/mx21;
vz=cross(x21,x22);mvz=sqrt(sum(vz.*vz));
ez=vz/mvz;
ey=cross(ez,ex);
T=[ex';ey';ez']';% transfprmation matrix

arwhl=0.20*mx21;% arrow head length
arwhang=15*pi/180;% arrow head angle, degree

P1=[0;0;0];
r=arwhl*sin(arwhang);
theta=linspace(0,2*pi,20);%resolution  can be increased
Pb=[
    arwhl*ones(size(theta))
    r*cos(theta)
    r*sin(theta)
    ];


P1=T*P1+p2;Pb=T*Pb+p2;

plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],lineColor,'linewidth',.75),hold on
fill3(Pb(1,:),Pb(2,:),Pb(3,:),headColor,'linewidth',0.5)
for i=1:length(theta)-1
%         fill3([P1(1) Pb(1,i) Pb(1,i+1) P1(1)],[P1(2) Pb(2,i) Pb(2,i+1) P1(2)],...
%         [P1(3) Pb(3,i) Pb(3,i+1) P1(3)],headColor,'linewidth',0.5)
    fill3([P1(1) Pb(1,i) Pb(1,i+1) P1(1)],[P1(2) Pb(2,i) Pb(2,i+1) P1(2)],...
        [P1(3) Pb(3,i) Pb(3,i+1) P1(3)],headColor,'EdgeColor','none');%,'linewidth',2)
end