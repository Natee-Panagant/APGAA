function panel = mesh2panel(node,ele)
panel(:,:,1)=[node(ele(:,1),1) node(ele(:,2),1) node(ele(:,3),1) node(ele(:,4),1) node(ele(:,1),1)];
panel(:,:,2)=[node(ele(:,1),2) node(ele(:,2),2) node(ele(:,3),2) node(ele(:,4),2) node(ele(:,1),2)];
panel(:,:,3)=[node(ele(:,1),3) node(ele(:,2),3) node(ele(:,3),3) node(ele(:,4),3) node(ele(:,1),3)];

% panel(:,:,1)=[node(ele(:,2),1) node(ele(:,1),1) node(ele(:,4),1) node(ele(:,3),1) node(ele(:,2),1)];
% panel(:,:,2)=[node(ele(:,2),2) node(ele(:,1),2) node(ele(:,4),2) node(ele(:,3),2) node(ele(:,2),2)];
% panel(:,:,3)=[node(ele(:,2),3) node(ele(:,1),3) node(ele(:,4),3) node(ele(:,3),3) node(ele(:,2),3)];

% npanel = size(ele,1);
% panel = zeros(npanel,5,3);
% for i = 1:npanel
%     % Detect Left/Right side panels
%     dy = node(ele(i,2),2)-node(ele(i,1),2);
%     dz = node(ele(i,2),3)-node(ele(i,1),3);
%     if dy > 0 %Right side panels
%         panel(i,:,1)=[node(ele(i,1),1) node(ele(i,2),1) node(ele(i,3),1) node(ele(i,4),1) node(ele(i,1),1)];
%         panel(i,:,2)=[node(ele(i,1),2) node(ele(i,2),2) node(ele(i,3),2) node(ele(i,4),2) node(ele(i,1),2)];
%         panel(i,:,3)=[node(ele(i,1),3) node(ele(i,2),3) node(ele(i,3),3) node(ele(i,4),3) node(ele(i,1),3)];
%     elseif dy < 0 %Left side panels
%         panel(i,:,1)=[node(ele(i,2),1) node(ele(i,1),1) node(ele(i,4),1) node(ele(i,3),1) node(ele(i,2),1)];
%         panel(i,:,2)=[node(ele(i,2),2) node(ele(i,1),2) node(ele(i,4),2) node(ele(i,3),2) node(ele(i,2),2)];
%         panel(i,:,3)=[node(ele(i,2),3) node(ele(i,1),3) node(ele(i,4),3) node(ele(i,3),3) node(ele(i,2),3)];
%     else %Vertical panels
%         if dz > 0
%         	panel(i,:,1)=[node(ele(i,1),1) node(ele(i,2),1) node(ele(i,3),1) node(ele(i,4),1) node(ele(i,1),1)];
%             panel(i,:,2)=[node(ele(i,1),2) node(ele(i,2),2) node(ele(i,3),2) node(ele(i,4),2) node(ele(i,1),2)];
%             panel(i,:,3)=[node(ele(i,1),3) node(ele(i,2),3) node(ele(i,3),3) node(ele(i,4),3) node(ele(i,1),3)];
%         elseif dz < 0
%             panel(i,:,1)=[node(ele(i,2),1) node(ele(i,1),1) node(ele(i,4),1) node(ele(i,3),1) node(ele(i,2),1)];
%             panel(i,:,2)=[node(ele(i,2),2) node(ele(i,1),2) node(ele(i,4),2) node(ele(i,3),2) node(ele(i,2),2)];
%             panel(i,:,3)=[node(ele(i,2),3) node(ele(i,1),3) node(ele(i,4),3) node(ele(i,3),3) node(ele(i,2),3)];
%         else
%             error('Mesh error, Please check!');
%         end
%     end
% end
end

