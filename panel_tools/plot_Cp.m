function gobj = plot_Cp(rst)
PanelDat = rst.PanelDat;

%% Generate panel plotting data
node = PanelDat.Nodes;
ele = PanelDat.WingPanel;
x = [node(ele(:,1),1) , node(ele(:,2),1) , node(ele(:,3),1) , node(ele(:,4),1) , node(ele(:,1),1)]';
y = [node(ele(:,1),2) , node(ele(:,2),2) , node(ele(:,3),2) , node(ele(:,4),2) , node(ele(:,1),2)]';
z = [node(ele(:,1),3) , node(ele(:,2),3) , node(ele(:,3),3) , node(ele(:,4),3) , node(ele(:,1),3)]';

%% Plot Pressure Coefficient (Cp)
cl = jet(7);
%
for mode = 1:2
    switch mode
        case 1 %VLM plot
            Cp = rst.Cp_VLM;
            figure;hold on;
            
            pat1 = fill3(x,y,z,Cp');
            colormap(gca,cl);
            colorbar;
            title({'VLM';'Cp'});
            box on
            grid on
            axis equal;
            set(gcf,'windowstate','maximize');
            view(90,90);

            ax1 = gca;
            ax = ax1;
            pat = pat1;

        case 2 %DLM plot
            
            Cp = rst.Cp_DLM;
            k = rst.Flight_Condition.k;
            Nk = numel(k);
            patT = [];
            axT = [];
            
            for i = 1:Nk
                %%
                figure;hold on;
                c = real(Cp{i})';
                pat1 = fill3(x,y,z,c);
                xlabel('Xaxis (m)')
                ylabel('Yaxis (m)')
                zlabel('Zaxis (m)')
                box on
                grid on
                colormap(gca,cl)
                colorbar(gca);
                title({'DLM';['Cp-Magnitude (@k = ' num2str(k(i)) ')']});
                axis equal;
                set(gcf,'windowstate','maximize');
                view(90,90);
                ax1 = gca;

                %%
                figure;hold on;
                c = imag(Cp{i})';
                pat2 = fill3(x,y,z,c);
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
                title({'DLM';['Cp-Phase (@k = ' num2str(k(i)) ')']});
                axis equal;
                set(gcf,'windowstate','maximize');
                view(90,90);
                ax2 = gca;
                patT = [patT;[{pat1},{pat2}]];
                axT = [axT; [{ax1},{ax2}]];
            end
            ax = axT;
            pat = patT;
    end

    gobj.pat = pat;
    gobj.ax = ax;
end