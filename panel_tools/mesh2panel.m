function panel = mesh2panel(node,ele)
% Convert mesh data (node+element) into panel data (combined into one array)
panel(:,:,1) = [node(ele(:,1),1) node(ele(:,2),1) node(ele(:,3),1) node(ele(:,4),1) node(ele(:,1),1)];
panel(:,:,2) = [node(ele(:,1),2) node(ele(:,2),2) node(ele(:,3),2) node(ele(:,4),2) node(ele(:,1),2)];
panel(:,:,3) = [node(ele(:,1),3) node(ele(:,2),3) node(ele(:,3),3) node(ele(:,4),3) node(ele(:,1),3)];
end

