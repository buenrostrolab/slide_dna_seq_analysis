% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of mouse cerebellum slide-DNA-seq data for Fig. 4

function[] = analyze_human_colon_cancer_4_dna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 4

expt = 'human_colon_cancer_4_dna';
sample = 'human_colon_cancer_4_dna_200114_13'; % Fig. 4

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

min_cov = 0;
pass_cov = sum(counts(:,1:end-1),2) > min_cov;

beads_display = beads(in_puck,:);
beads_filt = beads(in_puck & pass_cov,:);
counts_filt = counts(in_puck & pass_cov,:);
sel_beads = full(in_puck & pass_cov);

% save selected beads

dlmwrite(sprintf('processed/%s/%s.sel_beads.txt',expt,sample),sel_beads);

disp(sprintf('%s: Converted coordinates to microns and filtered beads',sec2time(toc)))

%% Smooth in xy space for PCA

counts_filt_smo_xy = zeros(size(counts_filt));

knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',50);

for i=1:length(knn_idx_xy); disp(i)
    counts_filt_smo_xy(i,:) = mean(counts_filt(knn_idx_xy(i,:),:),1);
end

%% Calculate chromosome z-scores

% set permutation parameters

num_perm = 100;
perm_chr_counts = zeros([size(counts,1),max(bins.chr_ind),num_perm],'uint8');
num_frags = size(sparse_frags,1);
frag_chr_bin = bins.chr_ind(sparse_frags.Var2);

% run permutations

for i=1:num_perm; disp(i)
    
    perm_bins = frag_chr_bin(randperm(num_frags));
    [C,ia,ic] = unique([sparse_frags{:,1} perm_bins],'rows');
    
    inds = sub2ind(size(perm_chr_counts),C(:,1),C(:,2),repmat(i,size(C,1),1));
    perm_chr_counts(inds) = accumarray(ic,1);
    perm_chr_counts(sub2ind(size(perm_chr_counts),C(:,1),C(:,2),repmat(i,size(C,1),1))) = accumarray(ic,1);
     
end

% calculate z-scores

chr_counts =  zeros([size(counts,1),max(bins.chr_ind)]);
[C,ia,ic] = unique([sparse_frags{:,1} frag_chr_bin],'rows');
chr_counts(sub2ind(size(chr_counts),C(:,1),C(:,2))) = accumarray(ic,1);
zscores = (chr_counts-mean(perm_chr_counts,3))./std(single(perm_chr_counts),[],3);

%% Smooth z-scores

chr_counts_filt = chr_counts(in_puck,:);
perm_chr_counts_filt = perm_chr_counts(in_puck,:,:);
zscores_filt = zscores(in_puck,:);

chr_counts_filt_smo = zeros(size(chr_counts_filt));
perm_chr_counts_filt_smo = zeros(size(perm_chr_counts_filt));

knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',50);

for i=1:length(knn_idx_xy); disp(i)
    chr_counts_filt_smo(i,:,:) = mean(chr_counts_filt(knn_idx_xy(i,:),:),1);
    perm_chr_counts_filt_smo(i,:,:) = mean(perm_chr_counts_filt(knn_idx_xy(i,:),:,:),1);
end

zscores_filt_smo = (chr_counts_filt_smo-mean(perm_chr_counts_filt_smo,3))./std(single(perm_chr_counts_filt_smo),[],3);

%% k-means clustering and annotation

[coeff,score,latent,tsquared,explained,mu] = pca(zscores_filt_smo(:,1:24)');

normal_cluster = 1;
subclone_1_cluster = 2;
subclone_2_cluster = 3;

rng(5)
clusters = kmeans(coeff(:,1:4),3);
clusters_labeled(clusters==1) = normal_cluster;
clusters_labeled(clusters==2) = subclone_1_cluster;
clusters_labeled(clusters==3) = subclone_2_cluster;
dlmwrite(sprintf('processed/%s/%s.clusters.txt',expt,sample),clusters_labeled)

