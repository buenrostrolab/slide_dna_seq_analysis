% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Plot watershed cluster outlines on existing figure

function[] = plot_ws_outlines(beads_filt,ws_cluster_labels)

for i=1:max(ws_cluster_labels)
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
end

end