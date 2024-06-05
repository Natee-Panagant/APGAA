function AcPanelGen01
addpath(genpath(pwd));
figure(4),clf,hold on
AC=feval('AcModel01');
nSurf=size(AC,2);
for i=1:nSurf
    iLabel1=strcmp({AC(i).Label}, 'Fuselage');%Fuselage 
    iLabel2=strcmp({AC(i).Label}, 'Wing');%Wing
    iLabel3=strcmp({AC(i).Label}, 'Horizontal tail');%Horizontal tail
    iLabel4=strcmp({AC(i).Label}, 'Vertical tail');%Vertical tail
    if iLabel1
        AcSurf=FuselageGridGen(AC(i));
        surf(AcSurf.X,AcSurf.Y,AcSurf.Z,'facecolor',[0.5 0.5 0.5],'facealpha',0.15,...
            'edgecolor',0.8*ones(1,3))
    end
    if iLabel2|iLabel3|iLabel4
%         AcPanel=fullWingPanelGen(AC(i),[0 1]);
%%%%%%%%%%%%%%%%%%%%
        if size(AC(i).SubSurf,2)==0
            AcSurf=fullWingPanelGen(AC(i),[0 1]);
        elseif size(AC(i).SubSurf,2)>0
            [iInterval,iSubSurf]=getWingSec(AC(i));
            for j=1:length(iSubSurf)
                if iSubSurf(j)==0
                    AcSurf=fullWingPanelGen(AC(i),[iInterval(j) iInterval(j+1)]);
                else
                    AcSurf=splitWingPanelGen(AC(i),iSubSurf(j));
                end
            end
        end
        %%%%%%%%%%%%%%%%
    end
    clear AcSurf
end
view(60,30)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
% axis equal
%%%%%%%%% sub-function %%%%%%%
function [iInterval,iSubSurf]=getWingSec(AC)
iInterval=[0 1];
for i=1:size(AC.SubSurf,2)
    iInterval=[iInterval AC.SubSurf(i).SpanData(1:2)];
end
iInterval=sort(unique(iInterval));
for i=1:(length(iInterval)-1)
    iSubSurf(i)=0;% no flap
    for j=1:size(AC.SubSurf,2)
        if AC.SubSurf(j).SpanData(1)==iInterval(i)&AC.SubSurf(j).SpanData(2)==iInterval(i+1)
            iSubSurf(i)=j;
            break
        end
    end
end
%%%%%%%%%% End of file %%%%%%%%%%
