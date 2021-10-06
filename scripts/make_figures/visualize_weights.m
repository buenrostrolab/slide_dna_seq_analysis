function[fig] = visualize_weights(profiles,control,bins,ylims,cmap)

sel = bins.chr ~= "chrY" & bins.chr ~= "chrM";

profiles = profiles(bins.chr ~= "chrY" & bins.chr ~= "chrM",:);
control = control(bins.chr ~= "chrY" & bins.chr ~= "chrM",:);
bins = bins(bins.chr ~= "chrY" & bins.chr ~= "chrM",:);

% detect num profiles and bin size

num_profiles = size(profiles,2);
bin_size = mode(bins.bin_len);
num_chrs = max(bins.chr_ind);

% detect hg19 or mm10

if sum(string(bins.chr) == 'chr22') > 1
    autosomal_bins = 1:max(find(bins.chr_ind == 22));
else
    autosomal_bins = 1:max(find(bins.chr_ind == 19));
end

% get chromosome lengths

chr_lens = zeros(num_chrs,1);
for i=1:num_chrs
   chr_lens(i) = max(bins.bin_end(bins.chr_ind == i));
end

% initialize figure

p = tight_subplot(1,1,[0.05 0.05],[0.2./1 0.2./1],[0.125 0.025]);
axes(p(1));

% loop through profiles

for i=1:num_profiles
    scatter(1:size(bins,1),profiles(:,i),25,mod(bins.chr_ind,2),'.'); hold on;
end
colormap(cmap)

% mark chromosomes

text_y = ylims(2) + ylims(2)/6;
tick_y = [ylims(2) - ylims(2)/50, ylims(2) + ylims(2)/50];

count = 0;

for chr=1:num_chrs

    count = count+ceil((chr_lens(chr)/bin_size));
    
    plot(repmat([count],2,1),tick_y,'Color','black','Marker','none'); hold on;
    if chr == 1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'chr1','HorizontalAlignment','center','FontSize',10); hold on;
    elseif chr == num_chrs
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'X','HorizontalAlignment','center','FontSize',10); hold on;
    elseif chr > 1 & chr <= num_chrs-1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,sprintf('%d',chr),'HorizontalAlignment','center','FontSize',10); hold on;
    end
    
end

% set axis labels and ticks

xlim([0 size(bins,1)]);
ylim([ylims(1) ylims(2)])

xlabel('Genomic position (1 Mb)','FontSize',12)
ylabel({'PC weight'},'FontSize',12);

xticks([])
%yticks([0:2:ylims(2)])
a = get(gca,'YTickLabel'); set(gca,'fontsize',12)

for i=1:ylims(2)-1
    plot([0 size(bins,1)],[i i],'Color',[.5 .5 .5],'Marker','none','LineStyle',':'); hold on;
end

% set figure size

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 6 1.5];  