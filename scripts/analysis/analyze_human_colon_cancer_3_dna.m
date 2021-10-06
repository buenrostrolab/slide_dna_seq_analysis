% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of mouse cerebellum slide-DNA-seq data for Fig. 3

function[] = analyze_human_colon_cancer_3_dna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 3

expt = 'human_colon_cancer_3_dna';
sample = 'human_colon_cancer_3_dna_191204_19'; % Fig. 3

% load sample data and genomic bins

sparse_counts = readtable(sprintf('data/%s/%s.sparse_counts_1Mb.txt',expt,sample));
counts = sparse(sparse_counts{:,1},sparse_counts{:,2},sparse_counts{:,3});
sparse_frags = readtable(sprintf('data/%s/%s.sparse_frags_1Mb.txt',expt,sample));
beads = readtable(sprintf('data/%s/%s.bead_locations.csv',expt,sample));
bins = readtable('reference/genomic_bins/hg19_1Mb_bins.txt');

disp(sprintf('%s: Loaded preprocessed data',sec2time(toc)))

%% Convert bead coordinates and filter stray beads

% convert to microns

beads.xcoord = beads.xcoord * (65/100);
beads.ycoord = beads.ycoord * (65/100);

% calculate weighted center and find beads more than 1.6 mm away

max_puck_radius = 1600;

weighted_center_x = mean(beads.xcoord); weighted_center_y = mean(beads.ycoord);
dist_from_center = pdist2([beads.xcoord,beads.ycoord],[weighted_center_x,weighted_center_y]);
in_puck = dist_from_center < max_puck_radius;

% filter beads with less than 100 reads

min_cov = 100;
pass_cov = sum(counts(:,1:end-1),2) > min_cov;

beads_display = beads(in_puck,:);
beads_filt = beads(in_puck & pass_cov,:);
counts_filt = counts(in_puck & pass_cov,:);
sel_beads = full(in_puck & pass_cov);

% save selected beads

dlmwrite(sprintf('processed/%s/%s.sel_beads.txt',expt,sample),sel_beads);

disp(sprintf('%s: Converted coordinates to microns and filtered beads',sec2time(toc)))

%%

[row col v] = find(counts_filt);
dlmwrite(sprintf('processed/%s/%s.sparse_counts_filt_1Mb.txt',expt,sample),[col row v], 'delimiter', ',')
writetable(beads_filt,sprintf('processed/%s/%s.beads_filt.txt',expt,sample))

%% Smooth in xy space for PCA

counts_filt_smo_xy = zeros(size(counts_filt));
filepath =  sprintf('processed/%s/%s.counts_filt_smo_xy.txt',expt,sample);

if isfile(filepath)
    counts_filt_smo_xy = dlmread(filepath);
else
    knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',50);
    for i=1:length(knn_idx_xy); disp(i)
        counts_filt_smo_xy(i,:) = mean(counts_filt(knn_idx_xy(i,:),:),1);
    end
    dlmwrite(filepath,counts_filt_smo_xy)
end

%% Normalize counts and select bins

% select bins

counts_filt_smo_xy_norm = counts_filt_smo_xy./mean(counts_filt_smo_xy,2);
cv = std(counts_filt_smo_xy_norm,[],1)./mean(counts_filt_smo_xy_norm,1);
pc_bins = mean(counts_filt_smo_xy_norm,1)>0 & mean(counts_filt_smo_xy_norm,1)<5 & cv<1;

% visualize mean and CV of selected bins

figure;
subplot(2,1,1);plot(mean(counts_filt_smo_xy_norm(:,pc_bins),1));xlabel('Genome bins');ylabel('Mean');xlim([0 length(cv(pc_bins))])
subplot(2,1,2);plot(cv(pc_bins));xlabel('Genome bins');ylabel('Coefficent of variation');xlim([0 length(cv(pc_bins))])

%% Run PCA on smoothed data and project beads

% run PCA

[coeff,score,latent,tsquared,explained,mu] = pca(counts_filt_smo_xy_norm(:,pc_bins),'Centered',true,'economy',true,'NumComponents',100);

% project beads

counts_filt_norm = counts_filt./mean(counts_filt,2);
pc_scores = full(calculatePCs(counts_filt_norm(:,pc_bins)',coeff(:,1:10),'center'));

%% Smooth pc_scores in xy and PC space

pc_scores_smo_xy = zeros(size(pc_scores));
pc_scores_smo_pc = zeros(size(pc_scores));
pc_scores_smo_both = zeros(size(pc_scores));

knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',50);
knn_idx_pc = knnsearch(pc_scores,pc_scores,'K',50);

for i=1:length(knn_idx_pc); disp(i)
    pc_scores_smo_pc(i,:) = mean(pc_scores(knn_idx_pc(i,:),:),1);
end

for i=1:length(knn_idx_xy); disp(i)
    pc_scores_smo_xy(i,:) = mean(pc_scores(knn_idx_xy(i,:),:),1);
    pc_scores_smo_both(i,:) = mean(pc_scores_smo_pc(knn_idx_xy(i,1:10),:),1);
end

dlmwrite(sprintf('processed/%s/%s.pc_scores_smo_both',expt,sample),pc_scores_smo_both)

%% Visualize projected PC scores

num_pcs = 4;

figure;
for pc=1:num_pcs

    subplot(5,num_pcs,(num_pcs*0)+pc); 
    scatter(beads_filt.xcoord,beads_filt.ycoord,20,pc_scores(:,pc),'.') ;
    title(sprintf('PC %d: %.02f%%',pc,explained(pc)));
    colormap(parula); caxis([prctile(pc_scores(:,pc),10) prctile(pc_scores(:,pc),90)])
    
    subplot(5,num_pcs,(num_pcs*1)+pc); 
    scatter(beads_filt.xcoord,beads_filt.ycoord,20,pc_scores_smo_xy(:,pc),'.') ;
    title(sprintf('PC %d: XY smoothing',pc));
    colormap(parula); caxis([prctile(pc_scores_smo_xy(:,pc),10) prctile(pc_scores_smo_xy(:,pc),90)])
    
    subplot(5,num_pcs,(num_pcs*2)+pc); 
    scatter(beads_filt.xcoord,beads_filt.ycoord,20,pc_scores_smo_pc(:,pc),'.') ;
    title(sprintf('PC %d: PC smoothing',pc));
    colormap(parula); caxis([prctile(pc_scores_smo_pc(:,pc),10) prctile(pc_scores_smo_pc(:,pc),90)])
    
    subplot(5,num_pcs,(num_pcs*3)+pc); 
    scatter(beads_filt.xcoord,beads_filt.ycoord,20,pc_scores_smo_both(:,pc),'.') ;
    title(sprintf('PC %d: joint smoothing',pc));
    colormap(parula); caxis([prctile(pc_scores_smo_both(:,pc),10) prctile(pc_scores_smo_both(:,pc),90)])
    
    subplot(5,num_pcs,(num_pcs*4)+pc);
    scatter(1:size(coeff,1),coeff(:,pc)',1,mod(bins.chr_ind(pc_bins,:),2)); xlim([1 size(coeff,1)]); colormap(jet)
    title(sprintf('PC %d: %.02f%%',pc,explained(pc)));
    xlabel('Genome bins');ylabel('PC loading')
    
end

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 10 8];  

%% k-means clustering and annotation

% k-means clustering

normal_cluster = 1;
subclone_1_cluster = 2;
subclone_2_cluster = 3;

rng(1)
clusters = kmeans(pc_scores_smo_both,3);
clusters_labeled(clusters==1) = normal_cluster;
clusters_labeled(clusters==3) = subclone_1_cluster;
clusters_labeled(clusters==2) = subclone_2_cluster;
dlmwrite(sprintf('processed/%s/%s.clusters.txt',expt,sample),clusters_labeled)

%% Fig. 3b - PC1, PC2, and clonal assignment

% PC1

sel_scores = pc_scores_smo_both(:,1);

fig_3b_pc1 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_scores,jet,[prctile(sel_scores,10) prctile(sel_scores,90)],0,4);
saveas(fig_3b_pc1,sprintf('figures/fig3/fig_3b_pc1.png'))

% PC2

sel_scores = -pc_scores_smo_both(:,2);

fig_3b_pc2 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_scores,jet,[prctile(sel_scores,10) prctile(sel_scores,90)],0,4);
saveas(fig_3b_pc2,sprintf('figures/fig3/fig_3b_pc2.png'))

% clonal assignment

cm = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880];
fig_3b_clusters = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,clusters_labeled,cm,'',0,4);
saveas(fig_3b_clusters,sprintf('figures/fig3/fig_3b_clusters.png'))

%% Calculate z-scores (run once per chromosome/region)

clusters_labeled = dlmread(sprintf('processed/%s/%s.clusters.txt',expt,sample));

%sel_bins = 1445:1545; % chr8q
%sel_bins = find(bins.chr_ind == 15); % chr15
%sel_bins = find(bins.chr_ind == 20); % chr20

%sel_bins = 122:250; % chr1q
%sel_bins = 1338:1398; % chr7q

num_frags = size(sparse_frags,1);
num_perms = 100;

sel_frags  = sparse_frags(sparse_frags.Var2 >= min(sel_bins) & sparse_frags.Var2 <= max(sel_bins),:);
perm_counts = zeros(size(counts,1),num_perms);

for i=1:num_perms; disp(i)
    
    perm_beads = sparse_frags.Var1(randperm(num_frags,size(sel_frags,1)));
    [C,ia,ic] = unique(perm_beads);
    perm_counts(sub2ind(size(perm_counts),C,repmat(i,size(C,1),1))) = accumarray(ic,1);
    
end

bin_counts = full(sum(counts_filt(:,sel_bins),2));
perm_counts_filt = perm_counts(in_puck & pass_cov,:);
bin_counts_smo = zeros(size(bin_counts));
perm_counts_smo = zeros(size(perm_counts_filt));

k = 10;
knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',k);
for i=1:length(knn_idx_xy); disp(i)
    bin_counts_smo(i) = sum(bin_counts(knn_idx_xy(i,1:k)),1);
    perm_counts_smo(i,:) = sum(perm_counts_filt(knn_idx_xy(i,1:k),:),1);
end

zscores = (bin_counts_smo-nanmean(perm_counts_smo,2))./std(perm_counts_smo,[],2);

norm_mean = nanmean(zscores(clusters_labeled == 1));
zscores_norm = zscores - norm_mean;

%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr8q',expt,sample),zscores_norm)
%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr15',expt,sample),zscores_norm)
%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr20',expt,sample),zscores_norm)

%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr1q',expt,sample),zscores_norm)
%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr7q',expt,sample),zscores_norm)

%% Fig. 3c - Genomic region z-scores for chr8q, chr15, and chr20

zscores_norm_chr6 = dlmread(sprintf('processed/%s/%s.zscores_norm_chr8q',expt,sample));
zscores_norm_chr15 = dlmread(sprintf('processed/%s/%s.zscores_norm_chr15',expt,sample));
zscores_norm_chr20 = dlmread(sprintf('processed/%s/%s.zscores_norm_chr20',expt,sample));

pvals_chr6 = -log10(normcdf(abs(zscores_norm_chr6).*-1).*2).*sign(zscores_norm_chr6);
pvals_chr15 = -log10(normcdf(abs(zscores_norm_chr15).*-1).*2).*sign(zscores_norm_chr15);
pvals_chr20 = -log10(normcdf(abs(zscores_norm_chr20).*-1).*2).*sign(zscores_norm_chr20);

fig_3c_chr8q = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pvals_chr6,redblue,[-10 10],0,4);
saveas(fig_3c_chr8q,sprintf('figures/fig3/fig_3c_chr8q.png'))

fig_3c_chr15 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pvals_chr15,redblue,[-5 5],0,4);
saveas(fig_3c_chr15,sprintf('figures/fig3/fig_3c_chr15.png'))

fig_3c_chr20 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pvals_chr20,redblue,[-5 5],0,4);
saveas(fig_3c_chr20,sprintf('figures/fig3/fig_3c_chr20.png'))

%% Find, smooth, and sort high-coverage beads

% select beads with high-coverage and high/low PC1

read_thresh = 2000;
pc_thresh = 5;
counts_high = counts_filt(sum(counts_filt,2)>read_thresh & abs(score(:,1))>pc_thresh,:);
clusters_high = clusters_labeled(sum(counts_filt,2)>read_thresh & abs(score(:,1))>pc_thresh);
[sort_clusters sort_idx] = sort(clusters_high);

% smooth over 10 Mb window

smooth_counts = zeros(size(counts_high));
window_size = 10;

for i=1:size(bins,1)
    nearby_bins = max(1,i-window_size/2):min(size(bins,1),i+window_size/2);
    same_chr =  nearby_bins(bins.chr_ind(nearby_bins) == bins.chr_ind(i));
    smooth_counts(:,i) = mean(counts_high(:,same_chr),2); 
end

counts_high_norm = normr(smooth_counts);
autosomal_bins = 1:max(find(bins.chr_ind == 22));

% sort within cluster by PC1

sel_score = score(sum(counts_filt,2)>read_thresh & abs(score(:,1))>pc_thresh,:);
sel_cluster_labels = clusters_labeled(sum(counts_filt,2)>read_thresh & abs(score(:,1))>pc_thresh);
sort_pc1_idx = zeros(size(sort_idx));

count = 1;
for i=1:3
    cluster_pc1 = sel_score(sel_cluster_labels==i,1);
    cluster_pc1_idx = find(sel_cluster_labels==i);
    [tmp sort_order] = sort(cluster_pc1,'ascend');
    sort_pc1_idx(count:count+size(cluster_pc1,1)-1) = cluster_pc1_idx(sort_order);
    count = count+size(cluster_pc1,1);
end

% save 

dlmwrite(sprintf('processed/%s/%s.high_cov_beads',expt,sample),counts_high_norm(sort_pc1_idx,autosomal_bins))

%% Fig. 3e - high-coverage bead profiles

fig_3e = figure;
visualize_beads(log2(counts_high_norm(sort_pc1_idx,autosomal_bins)./nanmedian(counts_high_norm(sort_pc1_idx,autosomal_bins),2)),bins(autosomal_bins,:),colorsJDB(0,0,'solar_extra_gray'),[-2 2]);
saveas(fig_3e,sprintf('figures/fig3/fig_3e.svg'))

%% Fig. S13b - Differential genomic regions (run 10x script first to get watershed/projection labels)

%zscores_norm = dlmread(sprintf('processed/%s/%s.zscores_norm_chr1q',expt,sample));
%zscores_norm = dlmread(sprintf('processed/%s/%s.zscores_norm_chr7q',expt,sample));
%zscores_norm = dlmread(sprintf('processed/%s/%s.zscores_norm_chr20',expt,sample));

colors = distinguishable_colors(101); colors(4,:) = [];
ws_cluster_labels = dlmread(sprintf('processed/%s/%s.ws_cluster_labels.txt',expt,sample));
proj_labels = uint16(dlmread(sprintf('processed/%s/%s.proj_label.txt',expt,sample)));
sel = ws_cluster_labels > 0;

p_val = -log10(normcdf(abs(zscores_norm).*-1).*2).*sign(zscores_norm);

fig_s13b = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord(~sel),beads_filt.ycoord(~sel),100,[.7 .7 .7],'','',0,4);
visualize_puck(beads_filt.xcoord(sel),beads_filt.ycoord(sel),100,p_val(sel),redblue,[-5 5],0,4);

% overlay watershed clusters

num_ws_clusters = max(ws_cluster_labels);
for i=1:num_ws_clusters
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color',colors(proj_labels(i),:),'LineWidth',5); hold on;
end

%saveas(fig_s13b,sprintf('figures/fig3/fig_s13b_chr1q.png'))
%saveas(fig_s13b,sprintf('figures/fig3/fig_s13b_chr7q.png'))
%saveas(fig_s13b,sprintf('figures/fig3/fig_s13b_chr20.png'))
