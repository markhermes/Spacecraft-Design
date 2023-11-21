function plot_nodes(trans)

z = trans(:,3);
i = find(abs(z)<0.0001);

node_set = trans(i,:);

node_1_indx = find(node_set(:,1) > 0);
node_1_indx = node_1_indx(1);

node_2_indx = find(node_set(:,1) < 0);
node_2_indx = node_2_indx(1);

plot3(node_set(node_1_indx,1),node_set(node_1_indx,2),...
    node_set(node_1_indx,3),'og',MarkerFaceColor='g',MarkerEdgeColor='k');
plot3(node_set(node_2_indx,1),node_set(node_2_indx,2),...
    node_set(node_2_indx,3),'og',MarkerFaceColor='g',MarkerEdgeColor='k');

end