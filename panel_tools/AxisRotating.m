function [X,Y,Z] = AxisRotating(X,Y,Z,Alpha,AxisVec,RefPt)

ex=AxisVec/sqrt(sum(AxisVec.*AxisVec));
ux=ex(1);uy=ex(2);uz=ex(3);
c=cos(Alpha);s=sin(Alpha);
T=[
    c+ux^2*(1-c),ux*uy*(1-c)-uz*s,ux*uz*(1-c)+uy*s
    ux*uy*(1-c)+uz*s,c+uy^2*(1-c),uy*uz*(1-c)-ux*s
    uz*ux*(1-c)-uy*s,uz*uy*(1-c)+ux*s,c+uz^2*(1-c)
    ];

for i=1:size(X,1)
    for j=1:size(X,2)
        RijPrime=RefPt+T*([X(i,j);Y(i,j);Z(i,j)]-RefPt);
        X(i,j)=RijPrime(1);
        Y(i,j)=RijPrime(2);
        Z(i,j)=RijPrime(3);
    end
end

