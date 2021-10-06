function[] = analyze_human_colon_cancer_3_scwgs()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% 0. Load genomic bins

bins = readtable('reference/genomic_bins/hg19_1Mb_bins.txt');

%% 1. Load 10x data

expt_10x = 'human_colon_cancer_3_scwgs';
sample1_10x = 'human_colon_cancer_3_scwgs_200227_01';
sample2_10x = 'human_colon_cancer_3_scwgs_210402_02';

% load sample data

counts1 = dlmread(sprintf('data/%s/%s.counts_1Mb.txt',expt_10x,sample1_10x));
counts2 = dlmread(sprintf('data/%s/%s.counts_1Mb.txt',expt_10x,sample2_10x));
counts_10x = [counts1; counts2];

counts_norm = counts_10x./mean(counts_10x,2);
sel_beads_10x = std(counts_norm,[],2)<1;
counts_filt_10x = counts_10x(sel_beads_10x,:);

dlmwrite(sprintf('processed/%s/sel_beads_10x.txt',expt_10x),sel_beads_10x)
clusters_10x = dlmread(sprintf('processed/%s/%s.clusters.txt',expt_10x,expt_10x));

%% 2. Load puck data

expt_puck = 'human_colon_cancer_3_dna';
sample_puck = 'human_colon_cancer_3_dna_191204_19'; % Fig. 3

% load sample data

sparse_counts = readtable(sprintf('data/%s/%s.sparse_counts_1Mb.txt',expt_puck,sample_puck));
counts_puck = sparse(sparse_counts{:,1},sparse_counts{:,2},sparse_counts{:,3});
beads = readtable(sprintf('data/%s/%s.bead_locations.csv',expt_puck,sample_puck));

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
pass_cov = sum(counts_puck(:,1:end-1),2) > min_cov;

beads_display = beads(in_puck,:);
beads_filt = beads(in_puck & pass_cov,:);
counts_filt_puck = counts_puck(in_puck & pass_cov,:);
sel_beads_puck = full(in_puck & pass_cov);

clusters_puck = dlmread(sprintf('processed/%s/%s.clusters.txt',expt_puck,sample_puck));

%% Load tracks for normalization and select bins

gc_table = readtable(sprintf('reference/gc_content/hg19_1Mb_gc.txt'));
map_table = readtable('reference/mappability/hg19_1Mb_map.txt');
repli_table = readtable('reference/rep_timing/hg19_1Mb_rep.txt');

% select appropriate column

gc_track = gc_table.Var8;
map_track = map_table.Var7;
rep_track = repli_table.Var7;

% set selection thresholds

gc_thresh = 0.35;
map_thresh = 0.7;
rep_thresh = 9000; % this was an arbitrary number chosen by me to represent missing data, normal range is around -4 to 4

% select bins

sel_bins = gc_track > gc_thresh & map_track > map_thresh & rep_track < rep_thresh;
sel_bins_auto = bins.chr_ind(sel_bins) <= 22;

gc_track_sel = gc_track(sel_bins);
map_track_sel = map_track(sel_bins);
rep_track_sel = rep_track(sel_bins);

%% 3. Norm 10x

counts_sel_10x = counts_filt_10x(:,sel_bins);
sel_beads = nansum(counts_sel_10x,2)>0;

[counts_gc_norm gc_fit gc_fig] = lowess_norm(counts_sel_10x,gc_track_sel,sel_bins_auto,sel_beads,1,"GC fit (1)");
[counts_map_norm map_fit map_fig] = lowess_norm(counts_gc_norm,map_track_sel,sel_bins_auto,sel_beads,1,"Mappability fit (2)");
[counts_rep_norm rep_fit rep_fig] = lowess_norm(counts_map_norm,rep_track_sel,sel_bins_auto,sel_beads,1,"Repli-seq fit only (3)");

counts_norm_10x = counts_rep_norm;
counts_mean_10x = counts_norm_10x./nanmean(counts_norm_10x,2);

%% 4. Norm puck

counts_sel_puck = counts_filt_puck(:,sel_bins);
sel_beads = nansum(counts_sel_puck,2)>0;

[counts_gc_norm gc_fit gc_fig] = lowess_norm(counts_sel_puck,gc_track_sel,sel_bins_auto,sel_beads,1,"GC fit (1)");
[counts_map_norm map_fit map_fig] = lowess_norm(counts_gc_norm,map_track_sel,sel_bins_auto,sel_beads,1,"Mappability fit (2)");
[counts_rep_norm rep_fit rep_fig] = lowess_norm(counts_map_norm,rep_track_sel,sel_bins_auto,sel_beads,1,"Repli-seq fit only (3)");

counts_norm_puck = counts_rep_norm;
counts_mean_puck = counts_norm_puck./nanmean(counts_norm_puck,2);

%% 5. Run PCA on 10x

[coeff,pc_scores_10x,latent,tsquared,explained,mu] = pca(counts_mean_10x,'Centered',false,'economy',true,'NumComponents',100);

%% 6. Project puck

pc_scores = full(calculatePCs(counts_mean_puck',coeff(:,1:10),''));

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

%% 7. Visualize PCs on puck

num_pcs = 6;

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
    scatter(1:size(coeff,1),coeff(:,pc)',1,mod(bins.chr_ind(sel_bins,:),2)); xlim([1 size(coeff,1)]); colormap(jet)
    title(sprintf('PC %d: %.02f%%',pc,explained(pc)));
    xlabel('Genome bins');ylabel('PC loading')
    
end

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 10 8];  

%% 8. Used reduced dimensionality for matching

num_pcs = 5;
counts_impute_puck = pc_scores(:,1:num_pcs)*coeff(:,1:num_pcs)';
counts_impute_10x = pc_scores_10x(:,1:num_pcs)*coeff(:,1:num_pcs)';

cluster_counts_impute_10x = zeros(max(clusters_10x),size(counts_norm_10x,2));

for i=1:max(clusters_10x)
   cluster_counts_impute_10x(i,:) = nanmean(counts_impute_10x(clusters_10x==i,:)); 
end

%% find thresh for PC1

figure; histogram(pc_scores_smo_both(:,1))
figure;
scatter(beads_filt.xcoord,beads_filt.ycoord,20,pc_scores_smo_both(:,1),'.');

pc1_thresh = 53;
figure;
scatter(beads_filt.xcoord,beads_filt.ycoord,20,pc_scores_smo_both(:,1)>pc1_thresh,'.');

%% Watershed clustering

%gauss_size = 25;
gauss_size = 50;
prc_thresh = 70;
hmin = 5;
area_min = 20000;

% select tumor beads

tumor_beads = pc_scores_smo_both(:,1)>pc1_thresh;
tumor_coords_int = uint16(beads_filt{tumor_beads,2:3});

% make tumor image

tumor_img = zeros(max(tumor_coords_int(:,1))+100,max(tumor_coords_int(:,2)+100));
tumor_img(sub2ind(size(tumor_img),tumor_coords_int(:,1),tumor_coords_int(:,2))) = 1;
tumor_img_smo = imgaussfilt(tumor_img,gauss_size);
bw = (tumor_img_smo>prctile(tumor_img_smo(:),prc_thresh));
%figure; imshow(bw',[])

% distance transform

dist_img = bwdist(~bw);
dist_img = -dist_img;
dist_img = imhmin(dist_img,hmin);
%figure; imshow(dist_img,[])

% watershed

ws_labels = watershed(dist_img);
ws_labels(~bw) = 0;
ws_rgb = label2rgb(ws_labels,'jet',[.5 .5 .5]);
%figure; imshow(ws_rgb)

% filter small areas

ws_labels_filt = zeros(size(ws_labels));
regions = regionprops(ws_labels,'Area');
count = 1;
for i=1:size(regions,1)
    if regions(i).Area > area_min
        ws_labels_filt(ws_labels==i) = count;
        count = count+1;
    end
end
ws_rgb_filt = label2rgb(ws_labels_filt,'jet',[.5 .5 .5]);
figure; imshow(ws_rgb_filt)

% transfer labels to beads

colors = distinguishable_colors(101); colors(4,:) = [];
ws_cluster_labels = zeros(size(beads_filt,1),1);
ws_cluster_labels = ws_labels_filt(sub2ind(size(tumor_img),uint16(beads_filt.xcoord),uint16(beads_filt.ycoord)));
num_ws_clusters = max(ws_cluster_labels);

% get cluster centers

for i=1:num_ws_clusters
    ws_cluster_centers(i,:) = mean(beads_filt{ws_cluster_labels==i,2:3});
    ws_cluster_counts(i,:) = nanmean(counts_norm_puck(ws_cluster_labels==i,:));
    ws_cluster_counts_impute(i,:) = nanmean(counts_impute_puck(ws_cluster_labels==i,:));
    ws_cluster_scores(i,:) = nanmean(pc_scores(ws_cluster_labels==i,:));
end

%% manually adjust clusters

cluster_counts_impute_10x_adj = zeros(6,2577);
cluster_counts_impute_10x_adj = cluster_counts_impute_10x(1:6,:);
cluster_counts_impute_10x_adj(2,:) = cluster_counts_impute_10x(8,:);
cluster_counts_impute_10x_adj(5,:) = cluster_counts_impute_10x(5,:) + cluster_counts_impute_10x(7,:);

%% cluster-cluster correlation (using imputed counts)

sel = ws_cluster_labels > 0;

cluster_corr = corr(cluster_counts_impute_10x_adj',ws_cluster_counts_impute');
[max_corr max_cluster] = max(cluster_corr,[],1);
sel_cluster = cluster_corr == max_corr;

bead_max_cluster = zeros(size(ws_cluster_labels));
bead_max_cluster(ws_cluster_labels>0) = max_cluster(ws_cluster_labels(ws_cluster_labels>0));

dlmwrite(sprintf('processed/%s/%s.ws_cluster_labels.txt',expt_puck,sample_puck),ws_cluster_labels)
dlmwrite(sprintf('processed/%s/%s.proj_label.txt',expt_puck,sample_puck),max_cluster)

%% Fig. 3f - scWGS coverage profiles

% remove doublet cluster

pc_scores_10x_no2 = pc_scores(clusters_10x ~= 2,:);
clusters_10x_no2 = clusters_10x(clusters_10x ~= 2);

% number clusters

cluster_labels_10x = zeros(size(clusters_10x_no2));
cluster_labels_10x(clusters_10x_no2 == 5 | clusters_10x_no2 == 7) = 1;
cluster_labels_10x(clusters_10x_no2 == 1) = 2;
cluster_labels_10x(clusters_10x_no2 == 6) = 3;
cluster_labels_10x(clusters_10x_no2 == 4) = 4;
cluster_labels_10x(clusters_10x_no2 == 3) = 5;
cluster_labels_10x(clusters_10x_no2 == 8) = 6;

% sort cells within cluster by PC1

[sort_clusters sort_idx] = sort(cluster_labels_10x);
sort_pc1_idx = zeros(size(sort_idx));

count = 1;
for i=1:max(cluster_labels_10x)
    cluster_pc1 = pc_scores_10x_no2(cluster_labels_10x==i,1);
    cluster_pc1_idx = find(cluster_labels_10x==i);
    [tmp sort_order] = sort(cluster_pc1,'ascend');
    sort_pc1_idx(count:count+size(cluster_pc1,1)-1) = cluster_pc1_idx(sort_order);
    count = count+size(cluster_pc1,1);
end

% norm

autosomal_bins = 1:max(find(bins.chr_ind == 22));
counts_filt_norm = counts_filt_10x./mean(counts_filt_10x,2);
counts_filt_norm_no2 = counts_filt_norm(clusters_10x ~= 2,:);

fig_3f = figure;
visualize_beads(log2(counts_filt_norm_no2(sort_pc1_idx,autosomal_bins)./nanmedian(counts_filt_norm_no2(sort_pc1_idx,autosomal_bins),2)),bins(autosomal_bins,:),colorsJDB(0,0,'solar_extra_gray'),[-2 2])

saveas(fig_3f,sprintf('figures/fig3/fig_3f.svg'))
saveas(fig_3f,sprintf('figures/fig3/fig_s13a.svg'))

%% Fig. 3g - projection of scWGS data

fig_3g = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4); hold on;
visualize_puck(beads_filt.xcoord(sel),beads_filt.ycoord(sel),100,colors(bead_max_cluster(sel),:),[],[],0,4)

for i=1:num_ws_clusters
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
end

saveas(fig_3g,sprintf('figures/fig3/fig_3g.png'))
saveas(fig_3f,sprintf('figures/fig3/fig_s13b_proj.png'))

%% Fig. S13 - scWGS projection

no_y = bins.chr_ind <= 23;
counts_norm_10x_filt = counts_norm_10x(clusters_10x ~= 2,:);

profiles = zeros(size(counts_10x,2),1);
profiles(sel_bins,1) = nanmean(counts_norm_10x_filt(cluster_labels_10x==6,:),1);
fig_s13c_10x_6 = figure; visualize_coverage(profiles(no_y,:),ones(sum(no_y),1),bins(no_y,:),[0 8],colors(2,:),2,"mode");
saveas(fig_s13c_10x_6,sprintf('figures/fig3/fig_s13c_10x_6.svg'))

profiles = zeros(size(counts_10x,2),1);
profiles(sel_bins,1) = nanmean(counts_norm_10x_filt(cluster_labels_10x==2,:),1);
fig_s13c_10x_2 = figure; visualize_coverage(profiles(no_y,:),ones(sum(no_y),1),bins(no_y,:),[0 8],colors(1,:),2,"mode");
saveas(fig_s13c_10x_2,sprintf('figures/fig3/fig_s13c_10x_2.svg'))

profiles = zeros(size(counts_10x,2),1);
profiles(sel_bins,1) = nanmean(counts_norm_10x_filt(cluster_labels_10x==3,:),1);
fig_s13c_10x_3 = figure; visualize_coverage(profiles(no_y,:),ones(sum(no_y),1),bins(no_y,:),[0 8],colors(6,:),2,"mode");
saveas(fig_s13c_10x_3,sprintf('figures/fig3/fig_s13c_10x_3.svg'))

profiles = zeros(size(counts_10x,2),1);
profiles(sel_bins,1) = nanmean(counts_norm_puck(bead_max_cluster==2,:),1);
fig_s13c_dna_2 = figure; visualize_coverage(profiles(no_y,:),ones(sum(no_y),1),bins(no_y,:),[0 8],colors(2,:),2,"mode");
saveas(fig_s13c_dna_2,sprintf('figures/fig3/fig_s13c_dna_2.svg'))

profiles = zeros(size(counts_10x,2),1);
profiles(sel_bins,1) = nanmean(counts_norm_puck(bead_max_cluster==1,:),1);
fig_s13c_dna_3 = figure; visualize_coverage(profiles(no_y,:),ones(sum(no_y),1),bins(no_y,:),[0 8],colors(1,:),2,"mode");
saveas(fig_s13c_dna_3,sprintf('figures/fig3/fig_s13c_dna_3.svg'))

profiles = zeros(size(counts_10x,2),1);
profiles(sel_bins,1) = nanmean(counts_norm_puck(bead_max_cluster==6,:),1);
fig_s13c_dna_6 = figure; visualize_coverage(profiles(no_y,:),ones(sum(no_y),1),bins(no_y,:),[0 8],colors(6,:),2,"mode");
saveas(fig_s13c_dna_6,sprintf('figures/fig3/fig_s13c_dna_6.svg'))
