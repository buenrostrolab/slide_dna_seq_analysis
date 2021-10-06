% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of mouse liver metastases slide-DNA-seq data for Fig. 1, S5, S7, S9

function[] = analyze_mouse_liver_met_1_dna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 1

expt = 'mouse_liver_met_1_dna';
%sample = 'mouse_liver_met_1_dna_191114_06';
sample = 'mouse_liver_met_1_dna_191114_05'; % Fig. S9

% load sample data and genomic bins

sparse_counts = readtable(sprintf('data/%s/%s.sparse_counts_1Mb.txt',expt,sample));
counts = sparse(sparse_counts{:,1},sparse_counts{:,2},sparse_counts{:,3});
beads = readtable(sprintf('data/%s/%s.bead_locations.csv',expt,sample));
bins = readtable('reference/genomic_bins/mm10_1Mb_bins.txt');

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

%% Replace with registered bead coordinates

beads_filt = readtable(sprintf('data/%s/%s.bead_locations_reg.csv',expt,sample));
beads_display = readtable(sprintf('data/%s/%s.bead_locations_disp.csv',expt,sample));

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
[counts_tn5_norm tn5_fit tn5_fig] = lowess_norm(counts_rep_norm,tn5_track_sel,sel_bins_auto,sel_beads,1,"Tn5 bias fit (4)");

counts_norm = counts_rep_norm;

%% Fig. S5c - Visualize slide-DNA-seq normalization

fig_s5c = figure;
subplot(2,4,1)
title(sprintf('GC %% (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 8000]); xlabel('GC %'); ylabel('Counts')

subplot(2,4,2)
title(sprintf('Mappability (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 8000]); xlabel('Mappability'); ylabel('Counts')

subplot(2,4,3)
title(sprintf('Replication timing (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 8000]); xlabel('Replication timing score'); ylabel('Counts')

subplot(2,4,4)
title(sprintf('Tn5 bias (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 8000]); xlabel('Tn5 bias'); ylabel('Counts')

subplot(2,4,5)
title(sprintf('GC %% (%.02f)',corr(nansum(counts_norm(sel_beads,sel_bins_auto),1)',gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,nansum(counts_norm(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('GC %'); ylabel('Norm counts')

subplot(2,4,6)
title(sprintf('Mappability (%.02f)',corr(nansum(counts_norm(sel_beads,sel_bins_auto),1)',map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,nansum(counts_norm(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Mappability'); ylabel('Norm counts')

subplot(2,4,7)
title(sprintf('Replication timing score (%.02f)',corr(nansum(counts_norm(sel_beads,sel_bins_auto),1)',rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,nansum(counts_norm(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Replication timing score'); ylabel('Norm counts')

subplot(2,4,8)
title(sprintf('Tn5 bias (%.02f)',corr(nansum(counts_norm(sel_beads,sel_bins_auto),1)',tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,nansum(counts_norm(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Tn5 bias'); ylabel('Norm counts')

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 12 6];  
saveas(fig_s5c,sprintf('figures/fig1/fig_s5c.svg'))

%% Normalize bulk cerebellum WGS

bulk_table = readtable('data/bulk_wgs/mouse_tumor_clone_wgs_190514_01.1Mb_coverage.txt');
bulk = bulk_table.Var4;

bulk_sel = bulk(sel_bins);

[bulk_gc_norm gc_fit gc_fig] = lowess_norm(bulk_sel',gc_track_sel,sel_bins_auto,ones(1,1),1,"GC fit (1)");
[bulk_map_norm map_fit map_fig] = lowess_norm(bulk_gc_norm,map_track_sel,sel_bins_auto,ones(1,1),1,"Mappability fit (2)");
[bulk_rep_norm rep_fit rep_fig] = lowess_norm(bulk_map_norm,rep_track_sel,sel_bins_auto,ones(1,1),1,"Repli-seq fit (3)");
[bulk_tn5_norm tn5_fit tn5_fig] = lowess_norm(bulk_rep_norm,tn5_track_sel,sel_bins_auto,ones(1,1),1,"Tn5 bias fit (4)");

bulk_norm = bulk_rep_norm';

%% Fig. S5d - Visualize bulk normalization

fig_s5d = figure;
subplot(2,4,1);
title(sprintf('GC %% (%.02f)',corr(bulk_sel(sel_bins_auto),gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,bulk_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 40000]); xlabel('GC %'); ylabel('Counts')

subplot(2,4,2)
title(sprintf('Mappability (%.02f)',corr(bulk_sel(sel_bins_auto),map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,bulk_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 40000]); xlabel('Mappability'); ylabel('Counts')

subplot(2,4,3)
title(sprintf('Replication timing (%.02f)',corr(bulk_sel(sel_bins_auto),rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,bulk_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 40000]);  xlabel('Replication timing score'); ylabel('Counts')

subplot(2,4,4)
title(sprintf('Tn5 bias (%.02f)',corr(bulk_sel(sel_bins_auto),tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,bulk_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 40000]); xlabel('Tn5 bias'); ylabel('Counts')

subplot(2,4,5)
title(sprintf('GC %% (%.02f)',corr(bulk_norm(sel_bins_auto),gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,bulk_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('GC %'); ylabel('Norm counts')

subplot(2,4,6)
title(sprintf('Mappability (%.02f)',corr(bulk_norm(sel_bins_auto),map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,bulk_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Mappability'); ylabel('Norm counts')

subplot(2,4,7)
title(sprintf('Replication timing score (%.02f)',corr(bulk_norm(sel_bins_auto),rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,bulk_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Replication timing score'); ylabel('Norm counts')

subplot(2,4,8)
title(sprintf('Tn5 bias (%.02f)',corr(bulk_norm(sel_bins_auto),tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,bulk_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Tn5 bias'); ylabel('Norm counts')

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 12 6];  
saveas(fig_s5d,sprintf('figures/fig1/fig_s5d.svg'))

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
tumor_cluster = 2;

rng(2);
clusters = kmeans(pc_scores_smo_both,2);
clusters_labeled(clusters==1) = normal_cluster;
clusters_labeled(clusters==2) = tumor_cluster;
eva = evalclusters(pc_scores_smo_both,'kmeans','CalinskiHarabasz','KList',[2:6])
dlmwrite(sprintf('processed/%s/%s.clusters.txt',expt,sample),clusters_labeled)

%% Fig. 1h, Fig. S9a, Fig. S7e - PC1 and tumor/normal bead assignments

% PC1

pc1_smo = pc_scores_smo_both(:,1);

fig_pc = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,pc1_smo,jet,[prctile(pc1_smo,10) prctile(pc1_smo,90)],0,4);
%saveas(fig_pc,sprintf('figures/fig1/fig_1h.png')) % for rep 1 (sample 06)
%saveas(fig_pc,sprintf('figures/fig1/fig_s9a_pc1_rep1.png')) % for rep 1 (sample 06)
%saveas(fig_pc,sprintf('figures/fig1/fig_s9a_pc1_rep2.png')) % for rep 2 (sample 05)

% tumor/normal bead assignments

cm = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
fig_clust = figure; hold on;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,clusters_labeled,cm,'',0,4);
%saveas(fig_clust,sprintf('figures/fig1/fig_s9a_clusters_rep1.png')) % for rep 1 (sample 06)
%saveas(fig_clust,sprintf('figures/fig1/fig_s7e.png')) % for rep 1 (sample 06)
%saveas(fig_clust,sprintf('figures/fig1/fig_s9a_clusters_rep2.png')) % for rep 2 (sample 05)

%% Fig. 1j - Genomic coverage by tumor/normal

profiles = zeros(size(counts_filt,2),2);
profiles(sel_bins,1) = nanmean(counts_norm(clusters_labeled == 1,:),1);
profiles(sel_bins,2) = nanmean(counts_norm(clusters_labeled == 2,:),1);

control_profile = zeros(size(counts_filt,2),1);
control_profile(sel_bins) = bulk_norm;

no_y = bins.chr_ind <= 20;

cm = [0 0.4470 0.7410; 0.6350 0.0780 0.1840];
fig_1j = figure;
visualize_coverage(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 6],cm,2,"mode");
saveas(fig_1j,sprintf('figures/fig1/fig_1j.svg'))

%% Fig. S9b - Normal/tumor genomic coverage by rep

control_profile = zeros(size(counts_filt,2),1);
control_profile(sel_bins) = bulk_norm;
no_y = bins.chr_ind <= 20;

% normal

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm(clusters_labeled == 1,:),1);

cm = [0 0.4470 0.7410];
fig_s9b_norm = figure;
visualize_coverage(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 6],cm,2,"mode");
%saveas(fig_s9b_norm,sprintf('figures/fig1/fig_s9b_norm_rep1.svg')) % for rep 1 (sample 06)
%saveas(fig_s9b_norm,sprintf('figures/fig1/fig_s9b_norm_rep2.svg')) % for rep 2 (sample 05)

% tumor

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm(clusters_labeled == 2,:),1);

cm = [0.6350 0.0780 0.1840];

fig_s9b_tumor = figure;
visualize_coverage(profiles(no_y,:),control_profile(no_y),bins(no_y,:),[0 6],cm,2,"mode");
%saveas(fig_s9b_tumor,sprintf('figures/fig1/fig_s9b_tumor_rep1.svg')) % for rep 1 (sample 06)
%saveas(fig_s9b_tumor,sprintf('figures/fig1/fig_s9b_tumor_rep2.svg')) % for rep 2 (sample 05)

%% Fig. S7 - clonal analysis workflow

% Fig. S7a - % variance explained

num_pcs = 8;

fig_s7a = figure; 
plot(1:num_pcs,explained(1:num_pcs)); 
 
xticks([1:1:num_pcs]); yticks([0:1:3]); 
a = get(gca,'YTickLabel'); set(gca,'fontsize',12)
xlabel('PCs','FontSize',12); ylabel('% variance explained','FontSize',12);

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 1.75 1.75];  

saveas(fig_s7a,sprintf('figures/fig1/fig_s7a.svg'))

% Fig. S7b - PC weights

pc = 1; % adjust PC here

profiles = zeros(size(counts_filt,2),1);
profiles(pc_bins,1) = coeff(:,pc);

fig_s7b = figure;
visualize_weights(profiles,ones(size(counts,2),1),bins,[prctile(coeff(:,pc),0) prctile(coeff(:,pc),100)],jet);
saveas(fig_s7b,sprintf('figures/fig1/fig_s7b_pc%d.svg',pc))

% Fig. S7c left - PC scores raw

sel_score = pc_scores(:,pc);

fig_s7c_raw = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_score,jet,[prctile(sel_score,10) prctile(sel_score,90)],0,4);
saveas(fig_s7c_raw,sprintf('figures/fig1/fig_s7c_raw_pc%d.png',pc))

% Fig. S7c right - PC scores raw

sel_score = pc_scores_smo_both(:,pc);

fig_s7c_smo = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],'','',1,4);
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sel_score,jet,[prctile(sel_score,10) prctile(sel_score,90)],0,4);
saveas(fig_s7c_smo,sprintf('figures/fig1/fig_s7c_smo_pc%d.png',pc))

% Fig. S7d - cluster evaluation

fig_s7d = figure; 
plot(2:6,eva.CriterionValues/10000); 
 
xlim([2 6]); xticks([2:1:6]); yticks([0.6:0.2:1.2])
a = get(gca,'YTickLabel'); set(gca,'fontsize',12)
xlabel('k','FontSize',12); ylabel('Calinski-Harabasz values','FontSize',12);

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 2 1.75];
saveas(fig_s7d,sprintf('figures/fig1/fig_s7d.svg'))
