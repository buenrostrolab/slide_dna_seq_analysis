function[fig out_profiles] = visualize_coverage_overlay_region(profiles,control,bins,ylims,cm,sel_bins)

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
    autosomal_bins = 1:max(find(bins.chr_ind == 19));
    build = "mm10";
end

% get chromosome lengths

chr_lens = zeros(num_chrs,1);
for i=1:num_chrs
   chr_lens(i) = max(bins.bin_end(bins.chr_ind == i));
end

% normalize coverage

norm_profiles = zeros(size(profiles));
norm_profiles = profiles./(control);
norm_profiles(isnan(norm_profiles)) = NaN;
norm_profiles(isinf(norm_profiles)) = NaN;
norm_profiles(norm_profiles == 0) = NaN;


%scaled_profiles = norm_profiles./median(norm_profiles(autosomal_bins,:),1)*2;
for i=1:size(norm_profiles,2)
   [f,xi] = ksdensity(norm_profiles(autosomal_bins,i)); [tmp idx] = max(f); bin_mode = xi(idx); 
   scaled_profiles(:,i) = norm_profiles(:,i)./bin_mode*1;
end

out_profiles = scaled_profiles(sel_bins,:);

figure('visible','off');
p = tight_subplot(1,1,[0 0],[0 0.],[0 0]);
axes(p(1));

sel_bins
size(sel_bins,1)
num_profiles

all_points = zeros(size(sel_bins,1).*num_profiles,1);
all_colors = zeros(size(sel_bins,1).*num_profiles,3);

count = 0;
for i=1:num_profiles
    all_points((i-1).*size(sel_bins,1)+1:(i).*size(sel_bins,1)) = scaled_profiles(sel_bins,i);
    all_colors((i-1).*size(sel_bins,1)+1:(i).*size(sel_bins,1),:) = repmat(cm(i,:),size(sel_bins,1),1);
end

all_bins = repmat(1:size(sel_bins,1),1,num_profiles);
rand_idx = randperm(size(sel_bins,1).*num_profiles);
scatter(all_bins(rand_idx),all_points(rand_idx)',500,all_colors(rand_idx,:),'.'); hold on;


xlim([0 size(sel_bins,1)]); ylim([ylims(1) ylims(2)])


yticks([0:2:6]); xticks([])
a = get(gca,'YTickLabel');  
set(gca,'fontsize',12)

for i=1:ylims(2)-1
    plot([0 size(bins,1)],[i i],'Color',[.5 .5 .5],'Marker','none','LineStyle',':'); hold on;
end

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 4*size(sel_bins,1)./75 4];  