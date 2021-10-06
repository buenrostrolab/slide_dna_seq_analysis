% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of mouse cerebellum slide-DNA-seq data for Fig. 2, S7, S12 

function[] = analyze_mouse_liver_met_2_dna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 2

expt = 'mouse_liver_met_2_dna';
sample = 'mouse_liver_met_2_dna_200114_10';

% load sample data and bins

sparse_counts = readtable(sprintf('data/%s/%s.sparse_counts_1Mb.txt',expt,sample));
counts = sparse(sparse_counts{:,1},sparse_counts{:,2},sparse_counts{:,3});
sparse_frags = readtable(sprintf('data/%s/%s.sparse_frags_1Mb.txt',expt,sample));
beads = readtable(sprintf('data/%s/%s.bead_locations.csv',expt,sample));
bins = readtable('reference/genomic_bins/mm10_1Mb_bins.txt');

disp(sprintf('%s: Loaded preprocessed data',sec2time(toc)))

%% Convert bead coordinates and filter stray beads

% convert to microns

beads.xcoord = -beads.xcoord * (65/100);
beads.ycoord = -beads.ycoord * (65/100);

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

%% Load tracks for normalization and select bins

gc_table = readtable(sprintf('reference/gc_content/mm10_1Mb_gc.txt'));
map_table = readtable('reference/mappability/mm10_1Mb_map.txt');
rep_table = readtable('reference/rep_timing/mm10_1Mb_rep.txt');
tn5_table = readtable(sprintf('reference/tn5_bias/mm10_1Mb_tn5.txt'));

% select appropriate column

gc_track = gc_table.Var8;
map_track = map_table.Var7;
rep_track = rep_table.Var7;
tn5_track = tn5_table.Var7;

% set selection thresholds

gc_thresh = 0.35;
map_thresh = 0.7;
rep_thresh = 9000; % this was an arbitrary number chosen to represent missing data, normal range is around -4 to 4

% select bins

sel_bins = gc_track > gc_thresh & map_track > map_thresh & rep_track < rep_thresh;
sel_bins_auto = bins.chr_ind(sel_bins) <= 19;

gc_track_sel = gc_track(sel_bins);
map_track_sel = map_track(sel_bins);
rep_track_sel = rep_track(sel_bins);
tn5_track_sel = tn5_track(sel_bins);

%% Normalize slide-DNA-seq

counts_sel = counts_filt(:,sel_bins);
sel_beads = nansum(counts_sel,2)>0;

[counts_gc_norm gc_fit gc_fig] = lowess_norm(counts_sel,gc_track_sel,sel_bins_auto,sel_beads,1,"GC fit (1)");
[counts_map_norm map_fit map_fig] = lowess_norm(counts_gc_norm,map_track_sel,sel_bins_auto,sel_beads,1,"Mappability fit (2)");
[counts_rep_norm rep_fit rep_fig] = lowess_norm(counts_map_norm,rep_track_sel,sel_bins_auto,sel_beads,1,"Repli-seq fit (3)");

counts_norm = counts_rep_norm;

%% Normalize bulk blood WGS

control_table = readtable('data/bulk_wgs/mouse_blood_wgs_190514_01.1Mb_coverage.txt');
control = control_table.Var4;

control_sel = control(sel_bins);

[control_gc_norm gc_fit gc_fig] = lowess_norm(control_sel',gc_track_sel,sel_bins_auto,ones(1,1),1,"GC fit (1)");
[control_map_norm map_fit map_fig] = lowess_norm(control_gc_norm,map_track_sel,sel_bins_auto,ones(1,1),1,"Mappability fit (2)");
[control_rep_norm rep_fit rep_fig] = lowess_norm(control_map_norm,rep_track_sel,sel_bins_auto,ones(1,1),1,"Repli-seq fit (3)");

bulk_norm = control_rep_norm';

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
clone_a_cluster = 2;
clone_b_cluster = 3;

rng(1)
clusters = kmeans(pc_scores_smo_both,3);
clusters_labeled(clusters==1) = normal_cluster;
clusters_labeled(clusters==2) = clone_a_cluster;
clusters_labeled(clusters==3) = clone_b_cluster;
eva = evalclusters(pc_scores_smo_both,'kmeans','CalinskiHarabasz','KList',[2:6])
dlmwrite(sprintf('processed/%s/%s.clusters.txt',expt,sample),clusters_labeled)

%% Fig. 2b - PC1, PC2, and clonal assignment

% PC1

sel_scores = pc_scores_smo_both(:,1);
fig_2b_pc1 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_scores,jet,[prctile(sel_scores,10) prctile(sel_scores,90)],0,4);
saveas(fig_2b_pc1,sprintf('figures/fig2/fig_2b_pc1.png'))

% PC2

sel_scores = pc_scores_smo_both(:,2);
fig_2b_pc2 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_scores,jet,[prctile(sel_scores,10) prctile(sel_scores,90)],0,4);
saveas(fig_2b_pc2,sprintf('figures/fig2/fig_2b_pc2.png'))

% clonal assignment

cm = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880];
fig_2b_clusters = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,clusters_labeled,cm,'',0,4);
saveas(fig_2b_clusters,sprintf('figures/fig2/fig_2b_clusters.png'))
saveas(fig_2b_clusters,sprintf('figures/fig2/fig_s7j.png'))

%% Calculate z-scores (run once per chromosome/region)

z_bins = 975:999; % chr6  126-150 Mb
%z_bins = find(bins.chr_ind == 15);
%z_bins = find(bins.chr_ind == 19);

num_frags = size(sparse_frags,1);
num_perms = 100;

sel_frags  = sparse_frags(sparse_frags.Var2 >= min(z_bins) & sparse_frags.Var2 <= max(z_bins),:);
perm_counts = zeros(size(counts,1),num_perms);

for i=1:num_perms; disp(i)
    
    perm_beads = sparse_frags.Var1(randperm(num_frags,size(sel_frags,1)));
    [C,ia,ic] = unique(perm_beads);
    perm_counts(sub2ind(size(perm_counts),C,repmat(i,size(C,1),1))) = accumarray(ic,1);
    
end

bin_counts = full(sum(counts_filt(:,z_bins),2));
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

dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr6',expt,sample),zscores_norm)
%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr15',expt,sample),zscores_norm)
%dlmwrite(sprintf('processed/%s/%s.zscores_norm_chr19',expt,sample),zscores_norm)

%% Fig. 2c - Genomic region z-scores for chr6 (126-150 Mb), chr15, and chr 19

zscores_norm_chr6 = dlmread(sprintf('processed/%s/%s.zscores_norm_chr6',expt,sample));
zscores_norm_chr15 = dlmread(sprintf('processed/%s/%s.zscores_norm_chr15',expt,sample));
zscores_norm_chr19 = dlmread(sprintf('processed/%s/%s.zscores_norm_chr19',expt,sample));

pvals_chr6 = -log10(normcdf(abs(zscores_norm_chr6).*-1).*2).*sign(zscores_norm_chr6);
pvals_chr15 = -log10(normcdf(abs(zscores_norm_chr15).*-1).*2).*sign(zscores_norm_chr15);
pvals_chr19 = -log10(normcdf(abs(zscores_norm_chr19).*-1).*2).*sign(zscores_norm_chr19);

fig_2c_chr6 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pvals_chr6,redblue,[-5 5],0,4);
saveas(fig_2c_chr6,sprintf('figures/fig2/fig_2c_chr6.png'))

fig_2c_chr15 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pvals_chr15,redblue,[-5 5],0,4);
saveas(fig_2c_chr15,sprintf('figures/fig2/fig_2c_chr15.png'))

fig_2c_chr19 = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pvals_chr19,redblue,[-5 5],0,4);
saveas(fig_2c_chr19,sprintf('figures/fig2/fig_2c_chr19.png'))

%% Fig. 2d - Genomic coverage by clone

profiles = zeros(size(counts_filt,2),3);
profiles(sel_bins,1) = nanmean(counts_norm(clusters_labeled == 1,:),1);
profiles(sel_bins,2) = nanmean(counts_norm(clusters_labeled == 2,:),1);
profiles(sel_bins,3) = nanmean(counts_norm(clusters_labeled == 3,:),1);

control_profile = zeros(size(counts_filt,2),1);
control_profile(sel_bins) = bulk_norm;

no_y = bins.chr_ind <= 20;
cm = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880];

fig_2d = figure;
visualize_coverage(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 3],cm,1,"mode");

rectangle('Position',[1 0 80-1 3]); hold on; % chr1 1-80 Mb
rectangle('Position',[880 0 905-880 3]); hold on; % chr6  30-55 Mb
rectangle('Position',[975 0 999-975 3]); hold on; % chr6  126-150 Mb
rectangle('Position',[min(find(bins.chr_ind==15)) 0 max(find(bins.chr_ind==15))-min(find(bins.chr_ind==15)) 3]) % chr15
rectangle('Position',[min(find(bins.chr_ind==19)) 0 max(find(bins.chr_ind==19))-min(find(bins.chr_ind==19)) 3]) % chr19

saveas(fig_2d,sprintf('figures/fig2/fig_2d.svg'))

%% Fig. 2d - zoom regions

zoom_bins = 1:80; % chr1 1-80 Mb
%zoom_bins = 880:905; % chr6 30-55 Mb
%zoom_bins = 975:999; % chr6 126-150 Mb
%zoom_bins = find(bins.chr_ind==15)'; % chr15
%zoom_bins = find(bins.chr_ind==19)'; % chr19

fig_2d_zoom = figure;
[fig zoom_profiles] = visualize_coverage_overlay_region(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 3],cm,zoom_bins');
ranksum(zoom_profiles(:,1),zoom_profiles(:,2))
ranksum(zoom_profiles(:,1),zoom_profiles(:,3))

saveas(fig,sprintf('figures/fig2/fig_2d_chr1.svg'))
%saveas(fig,sprintf('figures/fig2/fig_2d_chr6_1.svg'))
%saveas(fig,sprintf('figures/fig2/fig_2d_chr6_2.svg'))
%saveas(fig,sprintf('figures/fig2/fig_2d_chr15.svg'))
%saveas(fig,sprintf('figures/fig2/fig_2d_chr19.svg'))

%% Fig. S12 - Validation of ploidy and copy number of metastatic clones

no_y = bins.chr_ind <= 20;
cm = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880];

control_profile = zeros(size(counts_filt,2),1);
control_profile(sel_bins) = bulk_norm;

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm(clusters_labeled == 3,:),1);

fig_s12e = figure;
visualize_coverage(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 9],cm(3,:),3,"mode");

saveas(fig_s12e,sprintf('figures/fig2/fig_s12e.svg'))

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm(clusters_labeled == 2,:),1);

fig_s12f = figure;
visualize_coverage(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 6],cm(2,:),2,"mode");

saveas(fig_s12f,sprintf('figures/fig2/fig_s12f.svg'))

%% Fig. S7 - clonal analysis workflow

% Fig. S7f - % variance explained

num_pcs = 8;

fig_s7f = figure; 
plot(1:num_pcs,explained(1:num_pcs)); 
 
xticks([1:1:num_pcs]); yticks([0:1:3]); 
a = get(gca,'YTickLabel'); set(gca,'fontsize',12)
xlabel('PCs','FontSize',12); ylabel('% variance explained','FontSize',12);

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 1.75 1.75];  

saveas(fig_s7f,sprintf('figures/fig2/fig_s7f.svg'))

%% Fig. S7g - PC weights

pc = 2; % adjust PC here

profiles = zeros(size(counts_filt,2),1);
profiles(pc_bins,1) = coeff(:,pc);

fig_s7g = figure;
visualize_weights(profiles,ones(size(counts,2),1),bins,[prctile(coeff(:,pc),0) prctile(coeff(:,pc),100)],jet);
saveas(fig_s7g,sprintf('figures/fig2/fig_s7g_pc%d.svg',pc))

% Fig. S7h left - PC scores raw

sel_score = pc_scores(:,pc);

fig_s7h_raw = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_score,jet,[prctile(sel_score,10) prctile(sel_score,90)],0,4);
saveas(fig_s7h_raw,sprintf('figures/fig2/fig_s7h_raw_pc%d.png',pc))

% Fig. S7h right - PC scores raw

sel_score = pc_scores_smo_both(:,pc);

fig_s7h_smo = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_score,jet,[prctile(sel_score,10) prctile(sel_score,90)],0,4);
saveas(fig_s7h_smo,sprintf('figures/fig2/fig_s7h_smo_pc%d.png',pc))

% Fig. S7i - cluster evaluation

fig_s7d = figure; 
plot(2:6,eva.CriterionValues/10000); 
 
xlim([2 6]); xticks([2:1:6]); yticks([1:0.2:1.7])
a = get(gca,'YTickLabel'); set(gca,'fontsize',12)
xlabel('k','FontSize',12); ylabel('Calinski-Harabasz values','FontSize',12);

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 2 1.75];
saveas(fig_s7d,sprintf('figures/fig2/fig_s7i.svg'))
