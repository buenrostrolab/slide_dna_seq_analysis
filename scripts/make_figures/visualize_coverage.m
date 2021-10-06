function[fig scaled_profiles] = visualize_coverage(profiles,control,bins,ylims,cmap,ploidy,norm)

% detect num profiles and bin size

num_profiles = size(profiles,2);
bin_size = mode(bins.bin_len);
num_chrs = max(bins.chr_ind);

% detect hg19 or mm10

build = "";
if sum(string(bins.chr) == 'chr22') > 1
    autosomal_bins = 1:max(find(bins.chr_ind == 22));
    build = "hg19";
else
    autosomal_bins = bins.chr_ind <= 19;
    build = "mm10";
end

% get chromosome lengths


chr_lens = zeros(num_chrs,1);
for i=[unique(bins.chr_ind)'];
   chr_lens(i) = max(bins.bin_end(bins.chr_ind == i));
end

% normalize coverage

norm_profiles = zeros(size(profiles));
norm_profiles = profiles./(control);
norm_profiles(isnan(norm_profiles)) = NaN;
norm_profiles(isinf(norm_profiles)) = NaN;
norm_profiles(norm_profiles == 0) = NaN;

if norm == 'median'
    scaled_profiles = norm_profiles./nanmedian(norm_profiles(autosomal_bins,:),1)*ploidy;
elseif norm == 'mode'
    for i=1:size(norm_profiles,2)
       [f,xi] = ksdensity(norm_profiles(autosomal_bins,i)); [tmp idx] = max(f); bin_mode = xi(idx); 
       scaled_profiles(:,i) = norm_profiles(:,i)./bin_mode*ploidy;
    end
end

% initialize figure

p = tight_subplot(1,1,[0.05 0.05],[0.2./1 0.2./1],[0.125 0.025]);
axes(p(1));

% loop through profiles

%for i=1:num_profiles
%    scatter(1:size(bins,1),scaled_profiles(:,i),25,cmap(i,:),'.'); hold on;
%end

% randomize order of points

all_points = zeros(size(scaled_profiles,1).*num_profiles,1);
all_colors = zeros(size(scaled_profiles,1).*num_profiles,3);

for i=1:num_profiles
    all_points((i-1).*size(scaled_profiles,1)+1:(i).*size(scaled_profiles,1)) = scaled_profiles(:,i);
    all_colors((i-1).*size(scaled_profiles,1)+1:(i).*size(scaled_profiles,1),:) = repmat(cmap(i,:),size(scaled_profiles,1),1);
end

all_bins = repmat(1:size(bins,1),1,num_profiles);
rand_idx = randperm(size(scaled_profiles,1).*num_profiles);
scatter(all_bins(rand_idx),all_points(rand_idx)',25,all_colors(rand_idx,:),'.'); hold on;

% mark chromosomes

text_y = ylims(2) + ylims(2)/6;
tick_y = [ylims(2) - ylims(2)/50, ylims(2) + ylims(2)/50];

count = 0;

for chr=1:num_chrs

    count = count+ceil((chr_lens(chr)/bin_size));
    
    plot(repmat([count],2,1),tick_y,'Color','black','Marker','none'); hold on;
    if chr == 1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'chr1','HorizontalAlignment','center','FontSize',10); hold on;
    elseif (chr == 23 & build == "hg19") | (chr == 20 & build == "mm10")
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'X','HorizontalAlignment','center','FontSize',10); hold on;
    elseif (chr == 24 & build == "hg19") | (chr == 21 & build == "mm10")
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'Y','HorizontalAlignment','center','FontSize',10); hold on;
    elseif chr > 1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,sprintf('%d',chr),'HorizontalAlignment','center','FontSize',10); hold on;
    end
    
end

% set axis labels and ticks

xlim([0 size(bins,1)]);
ylim([ylims(1) ylims(2)])

xlabel('Genomic position (Mb)','FontSize',12)
ylabel({'Normalized','copy number'},'FontSize',12);

xticks([])
if ylims(2) ==  3 | ylims(2) ==  6
    yticks([0:ylims(2)./3:ylims(2)])
else
    yticks([0:2:ylims(2)])
end
a = get(gca,'YTickLabel'); set(gca,'fontsize',12)

for i=1:ylims(2)-1
    plot([0 size(bins,1)],[i i],'Color',[.5 .5 .5],'Marker','none','LineStyle',':'); hold on;
end

% set figure size

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 6 1.5];  