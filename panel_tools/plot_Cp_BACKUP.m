function gobj = plot_Cp(PanelDat,Cp,varargin)
if numel(varargin)>0
    k = varargin{1};
    mode = 2;
else
    mode = 1;
end

%% Generate panel plotting data
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
x = [node(ele(:,1),1) , node(ele(:,2),1) , node(ele(:,3),1) , node(ele(:,4),1) , node(ele(:,1),1)]';
y = [node(ele(:,1),2) , node(ele(:,2),2) , node(ele(:,3),2) , node(ele(:,4),2) , node(ele(:,1),2)]';
z = [node(ele(:,1),3) , node(ele(:,2),3) , node(ele(:,3),3) , node(ele(:,4),3) , node(ele(:,1),3)]';

%% Plot Pressure Coefficient (Cp)
cl = jet(7);
%
switch mode
    case 1 %VLM plot
        figure;hold on;
        %%
        pat1 = fill3(x,y,z,Cp');
        colormap(gca,cl);
        colorbar;
        title(['VLM - Cp']);
%         axis padded;
        axis equal;
%         view(-45,45);
        
        ax1 = gca;
        ax = ax1;
        pat = pat1;
        
    case 2 %DLM plot
        patT = [];
        axT = [];
        Nk = numel(Cp);
        for i = 1:Nk
            figure;hold on;
            
            %%
            subplot(1,2,1);
            c = real(Cp{i})';
            %             pat1 = fill3(x,y,z,real(Cp{i})');
            pat1 = fill3(x,y,z,c);
            xlabel('Xaxis (m)')
            ylabel('Yaxis (m)')
            zlabel('Zaxis (m)')
            box on
            grid on
            colormap(gca,cl)
            colorbar(gca);
            title(['DLM - Magnitude of Cp @k = ' num2str(k(i))]);
%             axis padded;
            axis equal;
%             view(-45,45);
            view(0,90)
            ax1 = gca;
            
            %%
            subplot(1,2,2);
            c = imag(Cp{i})';
            pat2 = fill3(x,y,z,c);
            %             pat2 = fill3(x,y,z,imag(Cp{i})');
            xlabel('Xaxis (m)')
            ylabel('Yaxis (m)')
            zlabel('Zaxis (m)')
            box on
            grid on
            colormap(gca,cl);
            colorbar(gca);
            climit = [min(c) max(c)];
            if abs(climit(1))>climit(2)
                colormap(gca,flipud(cl));
                colorbar(gca,'Direction','reverse');
            else
                colorbar(gca);
            end
            title(['DLM - Phase of Cp @k = ' num2str(k(i))]);
%             axis padded;
            axis equal;
%             view(-45,45);
            view(0,90)
            ax2 = gca;
            %
            patT = [patT;[{pat1},{pat2}]];
            axT = [axT; [{ax1},{ax2}]];
        end
        ax = axT;
        pat = patT;
        %         set(pat1,'Parent',hg);
        %         set(pat2,'Parent',hg);
        %
        %         set(hg,'Parent',ax);
        %         set(ax1,'Parent',hg);
        %         set(ax2,'Parent',hg);
        
end

gobj.pat = pat;
gobj.ax = ax;

