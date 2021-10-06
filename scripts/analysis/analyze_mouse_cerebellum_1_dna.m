% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of mouse cerebellum slide-DNA-seq data for Fig. 1, S3, S5, S6

function[] = analyze_mouse_cerebellum_1_dna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 1, S3, S5, S6, S10

expt = 'mouse_cerebellum_1_dna';
sample = 'mouse_cerebellum_1_dna_200114_14'; % Fig. 1, S3, S5, S6

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
sel_beads = dist_from_center < max_puck_radius;

% filter stray beads

beads_filt = beads(sel_beads,:);
counts_filt = counts(sel_beads,:);

% save selected beads

dlmwrite(sprintf('processed/%s/%s.sel_beads.txt',expt,sample),sel_beads);

disp(sprintf('%s: Converted coordinates to microns and filtered beads',sec2time(toc)))

%% Fig. 1D - Spatial % mtDNA distribution for mouse cerebellum

figure; 
fig_1f = visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sum(counts_filt(:,end),2)./sum(counts_filt(:,1:end),2)*100,winter,[0 25],1,4);
saveas(fig_1f,sprintf('figures/fig1/fig_1d.png'))

disp(sprintf('%s: Saved Fig. 1d ',sec2time(toc)))

%% Fig. S3b - Nuclear fragments/mtDNA fragments/% mtDNA per bead

figure;
fig_s3b_left = visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,log10(sum(counts_filt(:,1:end-1),2)),jet,[1 3],1,4);
saveas(fig_s3b_left,sprintf('figures/fig1/fig_s3b_left.png'))

figure;
fig_s3b_mid = visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,log10(sum(counts_filt(:,end),2)),jet,[0 2],1,4);
saveas(fig_s3b_mid,sprintf('figures/fig1/fig_s3b_mid.png'))

figure; 
fig_s3b_right = visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,sum(counts_filt(:,end),2)./sum(counts_filt(:,1:end),2)*100,winter,[0 25],1,4);
saveas(fig_s3b_right,sprintf('figures/fig1/fig_s3b_right.png'))

disp(sprintf('%s: Saved Fig. S3b ',sec2time(toc)))

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

%% Fig. S5a - Visualize slide-DNA-seq normalization

fig_s5a = figure;
subplot(2,4,1)
title(sprintf('GC %% (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2000]); xlabel('GC %'); ylabel('Counts')

subplot(2,4,2)
title(sprintf('Mappability (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2000]); xlabel('Mappability'); ylabel('Counts')

subplot(2,4,3)
title(sprintf('Replication timing (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2000]); xlabel('Replication timing score'); ylabel('Counts')

subplot(2,4,4)
title(sprintf('Tn5 bias (%.02f)',corr(nansum(counts_sel(sel_beads,sel_bins_auto),1)',tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,nansum(counts_sel(sel_beads,:),1),5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2000]); xlabel('Tn5 bias'); ylabel('Counts')

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
saveas(fig_s5a,sprintf('figures/fig1/fig_s5a.svg'))

%% Normalize bulk cerebellum WGS

control_table = readtable('data/bulk_wgs/mouse_cerebellum_wgs_210409_01.1Mb_coverage.txt');
control = control_table.Var4;

control_sel = control(sel_bins);

[control_gc_norm gc_fit gc_fig] = lowess_norm(control_sel',gc_track_sel,sel_bins_auto,ones(1,1),1,"GC fit (1)");
[control_map_norm map_fit map_fig] = lowess_norm(control_gc_norm,map_track_sel,sel_bins_auto,ones(1,1),1,"Mappability fit (2)");
[control_rep_norm rep_fit rep_fig] = lowess_norm(control_map_norm,rep_track_sel,sel_bins_auto,ones(1,1),1,"Repli-seq fit (3)");
[control_tn5_norm tn5_fit tn5_fig] = lowess_norm(control_rep_norm,tn5_track_sel,sel_bins_auto,ones(1,1),1,"Tn5 bias fit (4)");

control_norm = control_rep_norm';

%% Fig. S5b - Visualize bulk normalization

fig_s5b = figure;
subplot(2,4,1);
title(sprintf('GC %% (%.02f)',corr(control_sel(sel_bins_auto),gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,control_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 5000]); xlabel('GC %'); ylabel('Counts')

subplot(2,4,2)
title(sprintf('Mappability (%.02f)',corr(control_sel(sel_bins_auto),map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,control_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 5000]); xlabel('Mappability'); ylabel('Counts')

subplot(2,4,3)
title(sprintf('Replication timing (%.02f)',corr(control_sel(sel_bins_auto),rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,control_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 5000]);  xlabel('Replication timing score'); ylabel('Counts')

subplot(2,4,4)
title(sprintf('Tn5 bias (%.02f)',corr(control_sel(sel_bins_auto),tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,control_sel,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 5000]); xlabel('Tn5 bias'); ylabel('Counts')

subplot(2,4,5)
title(sprintf('GC %% (%.02f)',corr(control_norm(sel_bins_auto),gc_track_sel(sel_bins_auto)))); hold on;
scatter(gc_track_sel,control_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('GC %'); ylabel('Norm counts')

subplot(2,4,6)
title(sprintf('Mappability (%.02f)',corr(control_norm(sel_bins_auto),map_track_sel(sel_bins_auto)))); hold on;
scatter(map_track_sel,control_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Mappability'); ylabel('Norm counts')

subplot(2,4,7)
title(sprintf('Replication timing score (%.02f)',corr(control_norm(sel_bins_auto),rep_track_sel(sel_bins_auto)))); hold on;
scatter(rep_track_sel,control_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Replication timing score'); ylabel('Norm counts')

subplot(2,4,8)
title(sprintf('Tn5 bias (%.02f)',corr(control_norm(sel_bins_auto),tn5_track_sel(sel_bins_auto)))); hold on;
scatter(tn5_track_sel,control_norm,5,sel_bins_auto,'filled');
colormap(jet); ylim([0 2]); xlabel('Tn5 bias'); ylabel('Norm counts')

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 12 6];  
saveas(fig_s5b,sprintf('figures/fig1/fig_s5b.svg'))

%% Fig. 1f - Normalized genomic coverage for mouse cerebellum

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm,1);

control_profile = zeros(size(counts_filt,2),1);
control_profile(sel_bins,1) = control_norm;

auto_bins = bins.chr_ind < 20;

fig_1f = figure;
fig_1f = visualize_coverage(profiles(auto_bins,:),control_profile(auto_bins),bins(auto_bins,:),[0 6],[0 0.4470 0.7410],2,"median");
saveas(fig_1f,sprintf('figures/fig1/fig_1f.svg'))

disp(sprintf('%s: Saved Fig. 1f ',sec2time(toc)))

%% Fig. S6 - Normalization of copy number values

% Fig. S6a - slide-DNA-seq raw coverage

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_filt(:,sel_bins),1);

figure;
[fig_s6a_left scaled_profiles] = visualize_coverage(profiles(auto_bins,:),ones(sum(auto_bins),1),bins(auto_bins,:),[0 6],[0 0.4470 0.7410],2,"median");
title(sprintf('95 percentile range: %.2f - %.2f',prctile(scaled_profiles,2.5),prctile(scaled_profiles,97.5)))
saveas(fig_s6a_left,sprintf('figures/fig1/fig_s6a_left.svg'))

fig_s6a_right = figure; histogram(scaled_profiles(scaled_profiles>0),0:0.1:4)
xlim([0 4]); xticks([0:1:4]); ylim([0 700]); a = get(gca,'YTickLabel'); set(gca,'fontsize',10)
xlabel({'Copy number'},'FontSize',12); ylabel('# of bins','FontSize',12)
fig = gcf; fig.Units = 'inches'; fig.Position = [1 1 1.75 1.75];  
saveas(fig_s6a_right,sprintf('figures/fig1/fig_s6a_right.svg'))

% Fig. S6b - slide-DNA-seq norm coverage

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm,1);

figure;
[fig_s6b_left scaled_profiles] = visualize_coverage(profiles(auto_bins,:),ones(sum(auto_bins),1),bins(auto_bins,:),[0 6],[0 0.4470 0.7410],2,"median");
title(sprintf('95 percentile range: %.2f - %.2f',prctile(scaled_profiles,2.5),prctile(scaled_profiles,97.5)))
saveas(fig_s6b_left,sprintf('figures/fig1/fig_s6b_left.svg'))

fig_s6b_right = figure; histogram(scaled_profiles(scaled_profiles>0),0:0.1:4)
xlim([0 4]); xticks([0:1:4]); ylim([0 700]); a = get(gca,'YTickLabel'); set(gca,'fontsize',10)
xlabel({'Copy number'},'FontSize',12); ylabel('# of bins','FontSize',12)
fig = gcf; fig.Units = 'inches'; fig.Position = [1 1 1.75 1.75];  
saveas(fig_s6b_right,sprintf('figures/fig1/fig_s6b_right.svg'))

% Fig. S6c - bulk sequencing raw coverage

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = control(sel_bins);

figure;
[fig_s6c_left scaled_profiles] = visualize_coverage(profiles(auto_bins,:),ones(sum(auto_bins),1),bins(auto_bins,:),[0 6],[0 0.4470 0.7410],2,"median");
title(sprintf('95 percentile range: %.2f - %.2f',prctile(scaled_profiles,2.5),prctile(scaled_profiles,97.5)))
saveas(fig_s6c_left,sprintf('figures/fig1/fig_s6c_left.svg'))

fig_s6a_right = figure; histogram(scaled_profiles(scaled_profiles>0),0:0.1:4)
xlim([0 4]); xticks([0:1:4]); ylim([0 700]); a = get(gca,'YTickLabel'); set(gca,'fontsize',10)
xlabel({'Copy number'},'FontSize',12); ylabel('# of bins','FontSize',12)
fig = gcf; fig.Units = 'inches'; fig.Position = [1 1 1.75 1.75];  
saveas(fig_s6c_right,sprintf('figures/fig1/fig_s6c_right.svg'))

% Fig. S6d - bulk sequencing norm coverage

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = control_norm;

figure;
[fig_s6d_left scaled_profiles] = visualize_coverage(profiles(auto_bins,:),ones(sum(auto_bins),1),bins(auto_bins,:),[0 6],[0 0.4470 0.7410],2,"median");
title(sprintf('95 percentile range: %.2f - %.2f',prctile(scaled_profiles,2.5),prctile(scaled_profiles,97.5)))
saveas(fig_s6d_left,sprintf('figures/fig1/fig_s6d_left.svg'))

fig_s6d_right = figure; histogram(scaled_profiles(scaled_profiles>0),0:0.1:4)
xlim([0 4]); xticks([0:1:4]); ylim([0 700]); a = get(gca,'YTickLabel'); set(gca,'fontsize',10)
xlabel({'Copy number'},'FontSize',12); ylabel('# of bins','FontSize',12)
fig = gcf; fig.Units = 'inches'; fig.Position = [1 1 1.75 1.75];  
saveas(fig_s6d_right,sprintf('figures/fig1/fig_s6d_right.svg'))

% Fig. S6e - slide-DNA-seq / bulk norm coverage

profiles = zeros(size(counts_filt,2),1);
profiles(sel_bins,1) = nanmean(counts_norm,1);

figure;
[fig_s6e_left scaled_profiles] = visualize_coverage(profiles(auto_bins,:),control_profile(auto_bins),bins(auto_bins,:),[0 6],[0 0.4470 0.7410],2,"median");
title(sprintf('95 percentile range: %.2f - %.2f',prctile(scaled_profiles,2.5),prctile(scaled_profiles,97.5)))
saveas(fig_s6e_left,sprintf('figures/fig1/fig_s6e_left.svg'))

fig_s6e_right = figure; histogram(scaled_profiles(scaled_profiles>0),0:0.1:4)
xlim([0 4]); xticks([0:1:4]); ylim([0 700]); a = get(gca,'YTickLabel'); set(gca,'fontsize',10)
xlabel({'Copy number'},'FontSize',12); ylabel('# of bins','FontSize',12)
fig = gcf; fig.Units = 'inches'; fig.Position = [1 1 1.75 1.75];  
saveas(fig_s6e_right,sprintf('figures/fig1/fig_s6e_right.svg'))
