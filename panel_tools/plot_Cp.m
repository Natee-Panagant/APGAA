function plot_Cp(PanelDat,Cp,varargin)
if numel(varargin)>0
    k = varargin{1};
    mode = 2;
else
    mode = 1;
end

%Generate panel plotting data
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
x = [node(ele(:,1),1) , node(ele(:,2),1) , node(ele(:,3),1) , node(ele(:,4),1) , node(ele(:,1),1)]';
y = [node(ele(:,1),2) , node(ele(:,2),2) , node(ele(:,3),2) , node(ele(:,4),2) , node(ele(:,1),2)]';
z = [node(ele(:,1),3) , node(ele(:,2),3) , node(ele(:,3),3) , node(ele(:,4),3) , node(ele(:,1),3)]';

%Plot Pressure Coefficient (Cp)
switch mode
    case 1 %VLM plot
        figure;clf;hold on;
        fill3(x,y,z,Cp');
        colormap('jet');
        colorbar;
        title(['VLM - Cp']);
        axis equal;
        view(-30,30);
        set(gcf,'windowstate','maximize');
    case 2 %DLM plot
        Nk = numel(Cp);
        for i = 1:Nk
            figure;clf;hold on;

            subplot(1,2,1);
            fill3(x,y,z,real(Cp{i})');
            colormap('jet');
            colorbar;
            title(['DLM - Magnitude of Cp @k = ' num2str(k(1))]);
            axis equal;
            view(-30,30);
            set(gcf,'windowstate','maximize');

            subplot(1,2,2);
            fill3(x,y,z,imag(Cp{i})');
            colormap('jet');
            colorbar;
            title(['DLM - Phase of Cp @k = ' num2str(k(1))]);
            axis equal;
            view(-30,30);
            set(gcf,'windowstate','maximize');
        end
end
