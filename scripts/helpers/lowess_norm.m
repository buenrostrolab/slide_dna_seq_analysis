function[counts_norm fit_all fit_fig] = lowess_norm(counts,track,sel_fit,sel_beads,plot,plot_title)

    fit = smooth(track(sel_fit),full(nansum(counts(sel_beads,sel_fit),1))',0.1,'rlowess');
    fit_all = apply_smooth(fit,track(sel_fit),track);
    counts_norm = counts./fit_all';
    
    if plot == 1
        fit_fig = figure; 
        title(plot_title); hold on;
        scatter(track,nansum(counts(sel_beads,:),1),[],sel_fit,'filled');
        scatter(track,fit_all,'filled')
        colormap(jet); 
        ylim([0 prctile(nansum(counts(sel_beads,:),1),99)])
        
        fig = gcf;
        fig.Units = 'inches';
        fig.Position = [1 1 3 3];  
        
    else
       fit_fig = []; 
    end

end