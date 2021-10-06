% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of human colon cancer slide-RNA-seq data for Fig. 4, S15

function[] = analyze_human_colon_cancer_4_rna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 4

expt = 'human_colon_cancer_4_rna';
sample = 'human_colon_cancer_4_rna_200102_06';

% load sample data

sparse_counts = readtable(sprintf('data/%s/%s.sparse_expression.txt',expt,sample));
counts = sparse(sparse_counts{:,1},sparse_counts{:,2},sparse_counts{:,3})';
beads = readtable(sprintf('data/%s/%s.bead_locations.csv',expt,sample));
genes = readtable(sprintf('reference/gene_lists/human_genes.txt',expt,sample));

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

beads_filt = beads(in_puck,:);
counts_filt = counts(in_puck,:);

dlmwrite(sprintf('data/%s/%s.sel_beads.txt',expt,sample),full(in_puck));

%% Load bead cluster labels, UMAP coordinates, and registered beads

cluster_labels = ones(size(beads,1),1) * -1;
cluster_table = readtable(sprintf('processed/%s/%s.seurat_clusters.txt',expt,sample));
matched_beads = ismember(beads{:,1},cluster_table{:,1},'rows');
cluster_labels(matched_beads) = cluster_table{:,2};
cluster_labels_filt = cluster_labels(in_puck,:);

umap = ones(size(beads,1),2);
umap_table = readtable(sprintf('processed/%s/%s.seurat_umap.txt',expt,sample));
matched_beads = ismember(beads{:,1},umap_table{:,1},'rows');
umap(matched_beads,:) = umap_table{:,2:3};
umap_filt = umap(in_puck,:);

subclone_labels_filt = readtable(sprintf('processed/%s/%s.subclone_mask.txt',expt,sample));
subclone_labels_filt = subclone_labels_filt{:,3};

reg_bead_coords = dlmread(sprintf('processed/%s/%s.reg_beads.txt',expt,sample));

%% Filter beads without DNA assignment and few reads

read_thresh = 100;
dna_overlap = subclone_labels_filt>0 & sum(counts_filt,2)>read_thresh;

beads_display = beads_filt(subclone_labels_filt>0,:);

beads_filt = beads_filt(dna_overlap,:);
counts_filt = counts_filt(dna_overlap,:);

cluster_labels_filt = cluster_labels_filt(dna_overlap);
umap_filt = umap_filt(dna_overlap,:);
subclone_labels_filt = subclone_labels_filt(dna_overlap);
reg_bead_coords_filt = reg_bead_coords(dna_overlap,:);

writetable(beads_filt,sprintf('processed/%s/%s.beads_filt.txt',expt,sample))
writetable(beads_display,sprintf('processed/%s/%s.beads_display.txt',expt,sample))

%% Smooth in xy space 

filename = sprintf('processed/%s/%s.counts_filt_smo_xy.txt',expt,sample);

if isfile(filename)
    counts_filt_smo_xy = dlmread(filename);
else
    knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',50);
    tmp = full(counts_filt);

    counts_filt_smo_xy = zeros(size(counts_filt));
    for i = 1:length(knn_idx_xy);disp(i)
        counts_filt_smo_xy(i,:) = mean(tmp(knn_idx_xy(i,:),:),1);
    end

    dlmwrite(filename,counts_filt_smo_xy)
end

counts_filt_smo_xy_norm = normc(counts_filt_smo_xy')';

%% Label cell types

normal_beads = logical(sum(cluster_labels_filt == 1,2));
tumor_beads = logical(sum(cluster_labels_filt == [0 2 3 4],2));
immune_beads = logical(sum(cluster_labels_filt == [5:max(cluster_labels_filt)],2));

%normal_beads = logical(sum(cluster_labels_filt == [0 2 6],2));
%tumor_beads = logical(sum(cluster_labels_filt == [1 4],2));
%immune_beads = logical(sum(cluster_labels_filt == [3 5 7 8 9 10 11],2));

cell_type = zeros(size(normal_beads,1),1);
cell_type(normal_beads) = 1;
cell_type(tumor_beads) = 2;
cell_type(immune_beads) = 3;

dlmwrite(sprintf('processed/%s/%s.cell_type.txt',expt,sample),cell_type)

%% Watershed clustering

gauss_size = 50;
prc_thresh = 70;
hmin = 5;
area_min = 10000;

% select tumor beads

tumor_beads = logical(sum(cluster_labels_filt == [0 2 3 4],2));
tumor_coords_int = uint16(beads_filt{tumor_beads,2:3});

% make tumor image

tumor_img = zeros(max(tumor_coords_int(:,1))+5,max(tumor_coords_int(:,2)+5));
tumor_img(sub2ind(size(tumor_img),tumor_coords_int(:,1),tumor_coords_int(:,2))) = 1;
tumor_img_smo = imgaussfilt(tumor_img,gauss_size);
bw = (tumor_img_smo>prctile(tumor_img_smo(:),prc_thresh));
%figure; imshow(bw,[])

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
%figure; imshow(ws_rgb_filt)

% transfer labels to beads

colors = distinguishable_colors(101); colors(4,:) = [];
ws_cluster_labels = zeros(size(beads_filt,1),1);
ws_cluster_labels = ws_labels_filt(sub2ind(size(tumor_img),uint16(beads_filt.xcoord),uint16(beads_filt.ycoord)));
num_ws_clusters = max(ws_cluster_labels);

% get cluster centers

for i=1:num_ws_clusters
    ws_cluster_centers(i,:) = mean(beads_filt{ws_cluster_labels==i,2:3});
end

dlmwrite(sprintf('processed/%s/%s.watershed.txt',expt,sample),ws_cluster_labels)
dlmwrite(sprintf('processed/%s/%s.cluster_centers.txt',expt,sample),ws_cluster_centers)

%% Assign subclone labels to watershed clusters

ws_cluster_subclone_pct = zeros(num_ws_clusters,1);

for i=1:num_ws_clusters
    ws_cluster_subclone_pct(i) = sum(subclone_labels_filt(ws_cluster_labels==i) == 2) ./ sum(sum(subclone_labels_filt(ws_cluster_labels==i) == [2 3]));
end

subclone_thresh = 0.5;
ws_cluster_subclone_labels = ws_cluster_subclone_pct > subclone_thresh;
dlmwrite(sprintf('processed/%s/%s.subclone_labels.txt',expt,sample),ws_cluster_subclone_labels)

%% Calculate tumor and immune density

tumor_density = zeros(size(beads_filt,1),1);
immune_density = zeros(size(beads_filt,1),1);
normal_density = zeros(size(beads_filt,1),1);

knn_idx_xy = knnsearch(beads_filt{:,2:3},beads_filt{:,2:3},'K',50);

for i=1:length(knn_idx_xy); disp(i)
    tumor_density(i) = sum(tumor_beads(knn_idx_xy(i,:)))./50;
    immune_density(i) = sum(immune_beads(knn_idx_xy(i,:)))./50;
    normal_density(i) = sum(normal_beads(knn_idx_xy(i,:)))./50;
end

ws_cluster_tumor_neighbors = zeros(size(ws_cluster_subclone_labels));
ws_cluster_immune_neighbors = zeros(size(ws_cluster_subclone_labels));
ws_cluster_normal_neighbors = zeros(size(ws_cluster_subclone_labels));

for i=1:num_ws_clusters
   
    ws_cluster_tumor_neighbors(i) = mean(tumor_density(ws_cluster_labels==i));
    ws_cluster_immune_neighbors(i) = mean(immune_density(ws_cluster_labels==i));
    ws_cluster_normal_neighbors(i) = mean(normal_density(ws_cluster_labels==i));
    
end

dlmwrite(sprintf('processed/%s/%s.tumor_density.txt',expt,sample),tumor_density)

%% get counts for different cell types

ws_cluster_counts_all = zeros(num_ws_clusters,size(counts_filt,2));
ws_cluster_counts_tumor = zeros(num_ws_clusters,size(counts_filt,2));
ws_cluster_counts_immune = zeros(num_ws_clusters,size(counts_filt,2));

for i=1:num_ws_clusters; disp(i)
    ws_cluster_counts_all(i,:) = sum(counts_filt(ws_cluster_labels==i,:),1); 
end

sel_counts = ws_cluster_counts_all;
log_counts = log10(sel_counts+1);
norm_counts = quantilenorm(log_counts')';

%% Variance decomposition

sum_sq_mat = zeros(size(genes,1),4);
coeffs = zeros(size(genes,1),4);
pvals = zeros(size(genes,1),4);

% loop through genes

for i=1:size(genes,1); disp(i)

    % set x and y
    
    x = [ws_cluster_subclone_labels ws_cluster_tumor_neighbors ws_cluster_immune_neighbors];
    y = norm_counts(:,i);

    % fit stepwise and save coefficients
    
    p = stepwiselm(x,y,'linear','verbose',0,'Upper','linear');
    a = anova(p);
    in_model = p.VariableInfo.InModel; in_model(4) = 1;
    coeff_in_model = logical([1; in_model(1:3)]);
    sum_sq_mat(i,in_model) = a.SumSq;
    coeffs(i,coeff_in_model) = p.Coefficients.Estimate;
    pvals(i,in_model) = a.pValue;

end

dlmwrite(sprintf('processed/%s/%s.all_sum_sq_mat.txt',expt,sample),sum_sq_mat)
dlmwrite(sprintf('processed/%s/%s.all_coeffs.txt',expt,sample),coeffs)

%% Load preprocessed variance decomposition results

sum_sq_mat = dlmread(sprintf('processed/%s/%s.all_sum_sq_mat.txt',expt,sample));
coeffs = dlmread(sprintf('processed/%s/%s.all_coeffs.txt',expt,sample));

%% Process results

% convert to percent

prc_exp_mat = sum_sq_mat ./ sum(sum_sq_mat,2);

% combine density measurements

prc_exp_mat_comb(:,1) = prc_exp_mat(:,1);
prc_exp_mat_comb(:,2) = sum(prc_exp_mat(:,2:3),2);
prc_exp_mat_comb(:,3) = prc_exp_mat(:,4);

% add sign based on coefficients

prc_exp_mat_comb_dir = prc_exp_mat_comb;
prc_exp_mat_comb_dir(coeffs(:,2)<0,1) = prc_exp_mat_comb_dir(coeffs(:,2)<0,1) * -1;
prc_exp_mat_comb_dir(coeffs(:,3)-coeffs(:,4)<0,2) = prc_exp_mat_comb_dir(coeffs(:,3)-coeffs(:,4)<0,2) * -1;

dlmwrite(sprintf('processed/%s/%s.prc_exp_mat_comb_dir.txt',expt,sample),prc_exp_mat_comb_dir)

%% Fig. 4b - slide-RNA-seq tumor clusters

cm = [0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

fig_4b = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4); hold on;
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,cm(cell_type,:),[],[],0,4)

for i=1:num_ws_clusters
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
end

plot_label_backgrounds(ws_cluster_centers(:,1),ws_cluster_centers(:,2),125);
text(ws_cluster_centers(:,1),ws_cluster_centers(:,2),char(string(1:num_ws_clusters)'),'FontSize',14,'HorizontalAlignment','center')

saveas(fig_4b,sprintf('figures/fig4/fig_4b.png'))

%% Fig. 4c - slide-DNA-seq subclonal labels

fig_4c = figure;
cm = [0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(ws_cluster_labels>0),beads_filt.ycoord(ws_cluster_labels>0),100,ws_cluster_subclone_labels(ws_cluster_labels(ws_cluster_labels>0)),cm,[0 1],0,4)

for i=1:num_ws_clusters
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
end

saveas(fig_4c,sprintf('figures/fig4/fig_4c.png'))

%% Fig. 4d - slide-RNA-seq tumor cell density

fig_4d = figure;
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,tumor_density,colorsJDB(0,0,'brewer_spectra'),[0.5 1],0,4)

for i=1:num_ws_clusters
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
end

saveas(fig_4d,sprintf('figures/fig4/fig_4d.png'))

%% Fig. 4e - variance decomposition

sel = abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0;
den = datadensity(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),0.25);

fig_4e = figure;
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),10,den,'filled','MarkerEdgeColor',[.7 .7 .7]); hold on;
colormap((colorsJDB(0,0,'horizon')))

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

xlim([-1.1 1.1]); ylim([-1.1 1.1])
xticks([-1:1:1]); yticks([-1:1:1])
xticklabels({'100','0','100'})
yticklabels({'100','0','100'})
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;

sel_genes = ["LGALS3" "MCM7" "PLAG1" "PROM1"];
sel = logical(sum(string(genes.Var1)==sel_genes,2));
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),25,'r','filled');

text(prc_exp_mat_comb_dir(sel,1)*1.25,prc_exp_mat_comb_dir(sel,2)*1.25,char(sel_genes'),'HorizontalAlignment','center','FontSize',12)

ax = gca; ax.Position = [.15 .15 .75 .75];
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 4 4];

saveas(fig_4e,sprintf('figures/fig4/fig_4e.svg'))

%% Fig. 4f/g - selected gene exp

sel_beads = ws_cluster_labels>0;

fig_4f_mcm7 = figure;
exp = counts_filt_smo_xy_norm(:,string(genes.Var1) == "MCM7");
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],0,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_4f_mcm7,sprintf('figures/fig4/fig_4f_mcm7.png'))

fig_4f_plag1 = figure;
exp = counts_filt_smo_xy_norm(:,string(genes.Var1) == "PLAG1");
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],0,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_4f_plag1,sprintf('figures/fig4/fig_4f_plag1.png'))

fig_4g_lgals3 = figure;
exp = counts_filt_smo_xy_norm(:,string(genes.Var1) == "LGALS3");
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],0,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_4g_lgals3,sprintf('figures/fig4/fig_4g_lgals3.png'))

fig_4g_prom1 = figure;
exp = counts_filt_smo_xy_norm(:,string(genes.Var1) == "PROM1");
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],0,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_4g_prom1,sprintf('figures/fig4/fig_4g_prom1.png'))

%% Fig. 4i - variance decomposition of cell adhesion molecule binding

% plot all in gray

sel = abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0;

fig_4i = figure;
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),5,[.7 .7 .7],'filled','MarkerEdgeColor',[.7 .7 .7]); hold on;
colormap(colorsJDB(0,0,'horizon'))

ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
xlim([-1.1 1.1]); ylim([-1.1 1.1])
xticks([-1:1:1]); yticks([-1:1:1]) 
xticklabels({'100','0','100'});yticklabels({'100','0','100'})
ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;

% plot gene list in color and label selected genes

gene_list = readtable('reference/gene_lists/adhesion_gene_list.txt','ReadVariableNames',0);
sel_genes = ["CALD1","COL3A1","FLNB","ITGB2"];

sel = ismember(genes.Var1,gene_list.Var1) & (abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0);

den = datadensity(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),0.1);
[sort_den sort_idx] = sort(den);
sort_sel = find(sel); sort_sel = sort_sel(sort_idx);

scatter(prc_exp_mat_comb_dir(sort_sel,1),prc_exp_mat_comb_dir(sort_sel,2),10,sort_den,'filled','MarkerEdgeColor',[.7 .7 .7]);
sel = logical(sum(string(genes.Var1)==sel_genes,2));
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),25,'r','filled');
text(prc_exp_mat_comb_dir(sel,1)*1.25,prc_exp_mat_comb_dir(sel,2)*1.25,char(sel_genes'),'HorizontalAlignment','center','FontSize',10)

ax = gca; ax.Position = [0 0 1 1];
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 2 2];

saveas(fig_4i,sprintf('figures/fig4/fig_4i.svg'))
saveas(fig_4i,sprintf('figures/fig4/fig_s15h.svg'))

%% Fig. S15b - variance decomposition of E2F targets

% plot all in gray

sel = abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0;

fig_s15b = figure;
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),5,[.7 .7 .7],'filled','MarkerEdgeColor',[.7 .7 .7]); hold on;
colormap(colorsJDB(0,0,'horizon'))

ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
xlim([-1.1 1.1]); ylim([-1.1 1.1])
xticks([-1:1:1]); yticks([-1:1:1]) 
xticklabels({'100','0','100'});yticklabels({'100','0','100'})
ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;

% plot gene list in color and label selected genes

gene_list = readtable('reference/gene_lists/HALLMARK_E2F_TARGETS.txt','ReadVariableNames',0);
sel_genes = ["HMGB2","LYAR"];

sel = ismember(genes.Var1,gene_list.Var1) & (abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0);

den = datadensity(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),0.1);
[sort_den sort_idx] = sort(den);
sort_sel = find(sel); sort_sel = sort_sel(sort_idx);

scatter(prc_exp_mat_comb_dir(sort_sel,1),prc_exp_mat_comb_dir(sort_sel,2),10,sort_den,'filled','MarkerEdgeColor',[.7 .7 .7]);
sel = logical(sum(string(genes.Var1)==sel_genes,2));
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),25,'r','filled');
text(prc_exp_mat_comb_dir(sel,1)*1.25,prc_exp_mat_comb_dir(sel,2)*1.25,char(sel_genes'),'HorizontalAlignment','center','FontSize',10)

ax = gca; ax.Position = [0 0 1 1];
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 2 2];

saveas(fig_s15b,sprintf('figures/fig4/fig_s15b.svg'))

%% Fig. S15c - E2F target expression

gene_list = ["ANP32E","HMGB2","LMNB1","LYAR","MCM4","MCM7","NAA38","NAP1L1","RAN","TOP2A","USP1"];

fig_s15c = figure;
exp = sum(counts_filt_smo_xy_norm(:,ismember(string(genes.Var1),gene_list)),2);
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],1,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_s15c,sprintf('figures/fig4/fig_s15c.png'))

%% Fig. S15d - variance decomposition of MYC and MYC targets

% plot all in gray

sel = abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0;

fig_s15d = figure;
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),5,[.7 .7 .7],'filled','MarkerEdgeColor',[.7 .7 .7]); hold on;
colormap(colorsJDB(0,0,'horizon'))

ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
xlim([-1.1 1.1]); ylim([-1.1 1.1])
xticks([-1:1:1]); yticks([-1:1:1]) 
xticklabels({'100','0','100'});yticklabels({'100','0','100'})
ax.XAxis.FontSize = 8; ax.YAxis.FontSize = 8;

% plot gene list in color and label selected genes

gene_list = readtable('reference/gene_lists/HALLMARK_MYC_TARGETS_V1.txt','ReadVariableNames',0);
sel_genes = ["ILF2","MCM7","MYC","RRM1"];

sel = ismember(genes.Var1,gene_list.Var1) & (abs(prc_exp_mat_comb_dir(:,1)) > 0 | abs(prc_exp_mat_comb_dir(:,2)) > 0);

den = datadensity(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),0.1);
[sort_den sort_idx] = sort(den);
sort_sel = find(sel); sort_sel = sort_sel(sort_idx);

scatter(prc_exp_mat_comb_dir(sort_sel,1),prc_exp_mat_comb_dir(sort_sel,2),10,sort_den,'filled','MarkerEdgeColor',[.7 .7 .7]);
sel = logical(sum(string(genes.Var1)==sel_genes,2));
scatter(prc_exp_mat_comb_dir(sel,1),prc_exp_mat_comb_dir(sel,2),25,'r','filled');
text(prc_exp_mat_comb_dir(sel,1)*1.25,prc_exp_mat_comb_dir(sel,2)*1.25,char(sel_genes'),'HorizontalAlignment','center','FontSize',10)

ax = gca; ax.Position = [0 0 1 1];
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 2 2];

saveas(fig_s15d,sprintf('figures/fig4/fig_s15d.svg'))

%% Fig. S15e - MYC target expression

gene_list = ["CBX3";"CCNA2";"CCT3";"FBL";"HDGF";"HNRNPU";"ILF2";"MCM4";"MCM7";"NAP1L1";"PSMB3";"PTGES3";"RAN";"RRM1";"SLC25A3";"USP1"];

fig_s15e = figure;
exp = sum(counts_filt_smo_xy_norm(:,ismember(string(genes.Var1),gene_list)),2);
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],1,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_s15e,sprintf('figures/fig4/fig_s15e.png'))

gene_list = ["CBX3";"CCNA2";"CCT3";"FBL";"HDGF";"HNRNPU";"ILF2";"MCM4";"MCM7";"NAP1L1";"PSMB3";"PTGES3";"RAN";"RRM1";"SLC25A3";"USP1"];

fig_s15e_box = figure;
exp = sum(counts_filt_smo_xy_norm(:,ismember(string(genes.Var1),gene_list)),2);

sel_genes = ismember(string(genes.Var1),gene_list);
figure; notBoxPlot(sum(norm_counts(:,sel_genes),2),1-uint16(ws_cluster_subclone_labels))

xticklabels({'Subclone 1','Subclone 2'})
ylabel('Normalized spatial cluster expression','FontSize',12)
ylim([20 30]);
yticks([20:5:30])

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 3 3];

saveas(fig_s15e_box,sprintf('figures/fig4/fig_s15e_box.svg'))

%% Fig. S15f - MYC gene expression

fig_s15f = figure;
exp = counts_filt_smo_xy_norm(:,string(genes.Var1) == "MYC");
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],1,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_s15f,sprintf('figures/fig4/fig_15f.png'))

fig_s15f_scatter = figure;
scatter(ws_cluster_tumor_neighbors,norm_counts(:,string(genes.Var1) == 'MYC'),'filled'); hold on;

x = ws_cluster_tumor_neighbors;
y = norm_counts(:,string(genes.Var1) == 'MYC');

P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
plot(x,yfit,'r-.');

xlabel('Tumor cell density','FontSize',12)
ylabel('Normalized spatial cluster expression','FontSize',12)
ax.XAxis.FontSize = 10; ax.YAxis.FontSize = 10;

xlim([0.5 1])

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 3 3];
saveas(fig_s15f_scatter,sprintf('figures/fig4/fig_15f_scatter.svg'))

%% Fig. S15h - cell adhesion molecule binding expression

gene_list = ["AHNAK","CALD1","CHD11","CDH17","COL3A1","EPCAM","EVPL","FLNB","HDLBP","ITGB2","LASP1","PICALM","PTPN1","THY"];

fig_s15h = figure;
exp = sum(counts_filt_smo_xy_norm(:,ismember(string(genes.Var1),gene_list)),2);
visualize_puck(beads_display.xcoord,beads_display.ycoord,100,[.7 .7 .7],[],[],1,4)
visualize_puck(beads_filt.xcoord(sel_beads),beads_filt.ycoord(sel_beads),100,exp(sel_beads),colorsJDB(0,0,'brewer_spectra'),[prctile(exp(sel_beads),10) prctile(exp(sel_beads),90)],1,4);
plot_ws_outlines(beads_filt,ws_cluster_labels);
saveas(fig_s15h,sprintf('figures/fig4/fig_s15h.png'))

