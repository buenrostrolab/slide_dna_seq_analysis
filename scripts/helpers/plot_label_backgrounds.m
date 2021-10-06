function[] = plot_label_backgrounds(x,y,width)
    
    for i=1:size(x,1)
        rectangle('Position',[x(i)-(width/2) y(i)-(width/2) width width],...
            'Curvature',0.5,'FaceColor',[1 1 1],'EdgeColor','black','LineWidth',2); hold on;
    end

end
