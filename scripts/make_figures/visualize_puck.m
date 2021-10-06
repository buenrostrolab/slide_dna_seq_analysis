function[fig] = visualize_puck(x,y,size,col,cmap,col_lims,new_puck,img_size)

scatter(x,y,size,col,'.'); hold on;

weighted_center_x = nanmean(x); weighted_center_y = nanmean(y);
radius = 2000;

% set axes if new puck, otherwise overlay on previous

if new_puck == 1

    axis([(weighted_center_x - radius) (weighted_center_x + radius) (weighted_center_y - radius) (weighted_center_y + radius)]);
    set(gca,'xtick',[]); set(gca,'ytick',[]); axis off;
    
    scale_bar_x = weighted_center_x+1000; scale_bar_y = weighted_center_y-1500; 
    scale_bar_len = 500; scale_bar_height = 10;
    line([scale_bar_x scale_bar_x+scale_bar_len],[scale_bar_y scale_bar_y],'Color','black','LineWidth',3)

end

% set colormap

if ~isempty(cmap)
    colormap(cmap);
    if ~isempty(col_lims)
        caxis([min(col_lims) max(col_lims)]); 
    end
end

% set figure size

ax = gca; ax.Position = [0 0 1 1];
fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 img_size img_size];

end