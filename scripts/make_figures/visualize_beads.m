function[fig] = visualize_coverage(mat,bins,cmap,ylims)

% detect num profiles and bin size

bin_size = mode(bins.bin_len);
num_chrs = max(bins.chr_ind);

% detect hg19 or mm10

if num_chrs == 25
    autosomal_bins = 1:max(find(bins.chr_ind == 22));
elseif num_chrs == 22
    autosomal_bins = 1:max(find(bins.chr_ind == 19));
end

% get chromosome lengths

chr_lens = zeros(num_chrs,1);
for i=1:num_chrs
   chr_lens(i) = max(bins.bin_end(bins.chr_ind == i));
end

% plot beads

imagesc(mat); hold on;

% set colors

colormap(cmap)
caxis([ylims(1) ylims(2)])

% set axes

text_y = -size(mat,1)./25;
tick_y = 0;
a = get(gca,'YTickLabel');  
set(gca,'fontsize',12)
xticks([])
yticks([])

count = 1;
for chr=1:num_chrs; 

    count = count+ceil((chr_lens(chr)/bin_size));
    
    plot(repmat([count],2,1),[0 3000],'Color','black','Marker','none','LineWidth',2); hold on;
    if chr == 1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'chr1','HorizontalAlignment','center','FontSize',10); hold on;
    %elseif chr == num_chrs
    %    t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'22','HorizontalAlignment','center','FontSize',10); hold on;
    elseif chr >1 & mod(chr,2) == 1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,sprintf('%d',chr),'HorizontalAlignment','center','FontSize',10); hold on;
    %elseif chr == 9
    %    t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'...','HorizontalAlignment','center','FontSize',12); hold on;
    end

end

% set figure size

ax = gca; ax.Position = [.1 .1 .8 .8];
fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 4.5 3];  