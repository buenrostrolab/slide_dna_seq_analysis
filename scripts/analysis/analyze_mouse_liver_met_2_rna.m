% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Analysis of mouse liver metastases slide-RNA-seq data for Fig. 2, S11

function[] = analyze_mouse_liver_met_2_rna()
%% Initialize environment

tic
clearvars
home_dir = 'X:\\zchiang/slide_dna_seq2'; cd(home_dir);
addpath(genpath('scripts'))

disp(sprintf('%s: Initialized environment',sec2time(toc)))

%% Load data for Fig. 2

expt = 'mouse_liver_met_2_rna';
sample = 'mouse_liver_met_2_rna_201002_04';

% load sample data

sparse_counts = readtable(sprintf('data/%s/%s.sparse_expression.txt',expt,sample));
counts = sparse(sparse_counts{:,1},sparse_counts{:,2},sparse_counts{:,3})';
beads = readtable(sprintf('data/%s/%s.bead_locations.csv',expt,sample));
genes = readtable(sprintf('reference/gene_lists/mouse_genes.txt',expt,sample));

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

%% Load bead cluster labels and UMAP coordinates

cluster_labels = ones(size(beads,1),1) * -1;
cluster_table = readtable(sprintf('processed/%s/%s.seurat_clusters.txt',expt,sample),"ReadVariableNames",false);
matched_beads = ismember(beads{:,1},string(cluster_table{:,1}),'rows');
cluster_labels(matched_beads) = cluster_table{:,2};
cluster_labels_filt = cluster_labels(in_puck,:);

umap = ones(size(beads,1),2);
umap_table = readtable(sprintf('processed/%s/%s.seurat_umap.txt',expt,sample));
matched_beads = ismember(beads{:,1},umap_table{:,1},'rows');
umap(matched_beads,:) = umap_table{:,2:3};
umap_filt = umap(in_puck,:);

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

    dlmwrite(sprintf('processed/%s/%s/counts_filt_smo_xy.txt',expt,sample),counts_filt_smo_xy)
end

%% Watershed clustering

gauss_size = 50;
prc_thresh = 70;
hmin = 5;
area_min = 100000;

% select tumor beads

coords_int = uint16(beads_filt{:,2:3});
tumor_coords_int = uint16(beads_filt{tumor_beads,2:3});

% make tumor image

tumor_img = zeros(max(coords_int(:,1))+5,max(coords_int(:,2)+5));
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

ws_cluster_labels = zeros(size(beads_filt,1),1);
ws_cluster_labels = ws_labels_filt(sub2ind(size(tumor_img),uint16(beads_filt.xcoord),uint16(beads_filt.ycoord)));
dlmwrite(sprintf('processed/%s/%s.watershed.txt',expt,sample),ws_cluster_labels)

%% Fig. 2e

% load RCTD results and snRNA UMAP coords

rctd_table = readtable(sprintf('processed/%s/%s.rctd.csv',expt,sample));
umap = readtable(sprintf('processed/%s/%s.snrna_umap.csv',expt,sample));

% filter for singlets

beads_rctd = innerjoin(rctd_table,beads,'LeftKeys',1,'RightKeys',1);
beads_rctd_filt = beads_rctd(string(beads_rctd.spot_class) == 'singlet',:);

% manual annotation of clusters as normal, tumor, and immune

rctd_key = [14 3 7 11 10 13 5 9 6 4 12 1 2 8]';
type_key = [3 2 1 3 1 3 1 1 1 1 3 2 2 1];
type_cmap = [124 174 0; 248 118 109; 0 191 196]./255;

fig_2e = figure;
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,[.7 .7 .7],'','',1,8); hold on;
visualize_puck(beads_rctd_filt.xcoord,beads_rctd_filt.ycoord,100,type_cmap(type_key(ic),:),'','',0,4); 

% overlay watershed clusters

for i=1:max(ws_cluster_labels)
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
     ws_cluster_centers(i,:) = mean(beads_filt{ws_cluster_labels==i,2:3});
end

saveas(fig_2e,sprintf('figures/fig2/fig_2e.png'))

%% Perform differential expression analysis

counts_filt_full = full(counts_filt);

tumor_a = mean(counts_filt_full(ws_cluster_labels==1,:),1)';
tumor_b = mean(counts_filt_full(ws_cluster_labels==2,:),1)';

% permute cluster labels

num_perms = 100;
perm_exp = zeros(2,size(genes,1),num_perms);
in_cluster = find(ws_cluster_labels>0);
labels = ws_cluster_labels(in_cluster);
rng(1)

for i=1:num_perms; disp(i)

    perm_ws_cluster_labels(in_cluster) = labels(randperm(length(labels)));
    perm_exp(1,:,i) = mean(counts_full(perm_ws_cluster_labels==1,:),1)';
    perm_exp(2,:,i) = mean(counts_full(perm_ws_cluster_labels==2,:),1)';

end

% run z-test per gene

for i=1:size(genes,1); disp(i)
    [h,pvals(i),ci,zvals(i)] = ztest(log2(tumor_a(i)./tumor_b(i)),nanmean(log2(perm_exp(1,i,1:100)./perm_exp(2,i,1:100)),3)',std(log2(perm_exp(1,i,1:100)./perm_exp(2,i,1:100)),[],3)');
end

% multiple test correction

adj_pvals = pvals;
adj_pvals(isnan(pvals)) = 1;
adj_pvals = mafdr(adj_pvals,'BHFDR','true');
adj_pvals(adj_pvals==0) = normcdf(-abs(38.47)).*2;

dlmwrite(sprintf('processed/%s/%s.log2_exp.txt',expt,sample),-log2(nanmean(tumor_b,2)./nanmean(tumor_a,2)))
dlmwrite(sprintf('processed/%s/%s.adj_pvals.txt',expt,sample),adj_pvals)

%% Output differential gene expression table

out_table = table;
out_table.gene = genes.Var1;
out_table.clone1_exp = nanmean(tumor_a,2);
out_table.clone2_exp = nanmean(tumor_b,2);
out_table.log2_fc = -log2(nanmean(tumor_b,2)./nanmean(tumor_a,2));
out_table.fdr = adj_pvals';

sel_table = out_table(abs(out_table.log2_fc)>1 & (out_table.fdr <10^-2) & sum(counts_filt,1)'>100,:);
sort_table = sortrows(sel_table,5,'ascend');
writetable(sort_table,'tables/Supplementary_Table_4_mouse_clone_DEGs.txt');

%% Fig. 2f - volcano plot

log2_exp = dlmread(sprintf('processed/%s/%s.log2_exp.txt',expt,sample));
adj_pvals = dlmread(sprintf('processed/%s/%s.adj_pvals.txt',expt,sample));

fig_2f = figure;
scatter(log2_exp,-log10(adj_pvals),10,'b','filled','MarkerEdgeColor',[.7 .7 .7]); hold on; 
colormap(colorsJDB(0,0,'brewer_spectra'));

sel_genes = ["Hmga2" "Vim" "Tm4sf1" "Aqp5","S100a4","Vcan"];
sel = any(string(genes.Var1)==sel_genes,2);
scatter(log2_exp(sel),-log10(adj_pvals(sel)),25,'r','filled')
ordered = string(genes.Var1(any(string(genes.Var1)==sel_genes,2)));
text(log2_exp(sel),-log10(adj_pvals(sel)),ordered)

xlim([-8 8])
xticks(-8:4:8); yticks(0:100:300)
a = get(gca,'YTickLabel');  
set(gca,'fontsize',10)
xlabel('log2(B/A)','FontSize',12); 
ylabel('-log10(FDR)','FontSize',12)

ax = gca; ax.Position = [.15 .15 .75 .75];
fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 3 3];
saveas(fig_2f,sprintf('figures/fig2/fig_2f.svg'))

%% Fig. 2g - clone-specific gene expression

fig_2g_hmga2 = figure;
exp = counts_filt_smo_xy(:,string(genes.Var1) == "Hmga2");
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,exp,colorsJDB(0,0,'brewer_spectra'),[prctile(exp,10) prctile(exp,90)],1,4);
saveas(fig_2g_hmga2,sprintf('figures/fig2/fig_2g_hmga2.png'))
saveas(fig_2g_hmga2,sprintf('figures/fig2/fig_s11a_right.png'))

fig_2g_tm4sf1 = figure;
exp = counts_filt_smo_xy(:,string(genes.Var1) == "Tm4sf1");
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,exp,colorsJDB(0,0,'brewer_spectra'),[prctile(exp,10) prctile(exp,90)],1,4);
saveas(fig_2g_tm4sf1,sprintf('figures/fig2/fig_2g_tm4sf1.png'))

fig_2g_aqp5 = figure;
exp = counts_filt_smo_xy(:,string(genes.Var1) == "Aqp5");
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,exp,colorsJDB(0,0,'brewer_spectra'),[prctile(exp,10) prctile(exp,90)],1,4);
saveas(fig_2g_aqp5,sprintf('figures/fig2/fig_2g_aqp5.png'))

%% Fig. S11b - snRNA UMAP

% keep seurat coloring

rna_cmap = [248 118 109; 230 134 19; 205 150 0; 171 163 0; 124 174 0; 12 183 2; 0 190 103; 0 193 154; 0 191 196; 0 184 231; 0 169 255; 132 148 255; 199 124 255; 237 104 237; 255 97 204; 255 104 161]./255;
umap_key = [5 1 2 10 11 12 7 3 9 6 13 8 4 14]';

fig_s11b = figure; scatter(umap.UMAP_1,umap.UMAP_2,1,rna_cmap(umap_key(umap.Var4),:),'filled')
xlim([-15 15]); ylim([-15 15])

xlabel('UMAP1','FontSize',12); 
ylabel('UMAP2','FontSize',12)

ax = gca; ax.Position = [.2 .2 .7 .7];
fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 3 3];
saveas(fig_s11b,sprintf('figures/fig2/fig_s11b.svg'))

fig_s11b_leg = figure; 
for i=1:16
    scatter(1,1,100,rna_cmap(i,:),'filled'); hold on;
end
lgd = legend('Tumor I','Tumor II','Tumor III','Interferon','Hepatocyte I','Hepatocyte II','Hepatic Stellate Cell','VSMC','Ribosomal+',...
    'LSEC','Kupffer/Monocyte','Monocyte','T','B')
lgd.FontSize = 12;
saveas(fig_s11b_leg,sprintf('figures/fig2/fig_s11b_legend.svg'))

%% Fig. S11e - RCTD projection

% sort clusters by descending frequency

[clusters ia ic] = unique(beads_rctd_filt.first_type);
counts = accumarray(ic,1);
[sort_counts sort_order] = sort(counts,'descend');

% manual annotation of clusters as normal, tumor, and immune

rctd_key = [14 3 7 11 10 13 5 9 6 4 12 1 2 8]';
type_key = [3 2 1 3 1 3 1 1 1 1 3 2 2 1];
type_cmap = [124 174 0; 248 118 109; 0 191 196]./255;

fig_s11e = figure;
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,[.7 .7 .7],'','',1,4); hold on;
for i=[sort_order']
    sel = ic == i;
    visualize_puck(beads_rctd_filt.xcoord(sel),beads_rctd_filt.ycoord(sel),100,rna_cmap(rctd_key(ic(sel)),:),'','',0,4); 
end

% overlay watershed clusters

for i=1:max(ws_cluster_labels)
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
     ws_cluster_centers(i,:) = mean(beads_filt{ws_cluster_labels==i,2:3});
end

saveas(fig_s11e,sprintf('figures/fig2/fig_s11e.png'))

%% Fig. S11f - Monocyte localization

sel = string(beads_rctd_filt.first_type) == 'monocyte/DC';

fig_s11f = figure;
visualize_puck(beads_filt.xcoord,beads_filt.ycoord,100,[.7 .7 .7],'','',1,4); hold on;
visualize_puck(beads_rctd_filt.xcoord(sel),beads_rctd_filt.ycoord(sel),500,rna_cmap(12,:),'','',0,4); 

for i=1:max(ws_cluster_labels)
     x = beads_filt.xcoord(ws_cluster_labels==i);
     y = beads_filt.ycoord(ws_cluster_labels==i);
     k = boundary(x,y);
     plot(x(k),y(k),'Color','black','LineWidth',5); hold on;
     ws_cluster_centers(i,:) = mean(beads_filt{ws_cluster_labels==i,2:3});
end

saveas(fig_s11f,sprintf('figures/fig2/fig_s11f.png'))
