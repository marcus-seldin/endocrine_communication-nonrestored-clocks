#Required packages for analysis
library(dplyr)
library(reshape2)
library(viridis)
library(WGCNA) 
library(pheatmap)
library(ggplot2)
library(forcats)
library(colormap)
allowWGCNAThreads() 

#Change working directory as needed. All files should be in the same folder with this script.
setwd('path_to_files/')

#Load in population datasets and lists
#Read in curated lists of secreted proteins, along with the database of traits.
sec_prots = read.delim('secreted_proteins.txt')
adip = read.delim('HMDP_HFHS_trx_adipose.txt')
musc = read.delim('HMDP_HFHS_trx_muscle.txt')
liver = read.delim('HMDP_HFHS_trx_liver.txt')
int = read.delim('HMDP_HFHS_trx_intestine.txt')
hypo = read.delim('HMDP_HFHS_trx_hypothalamus.txt')
h_pathway = read.delim('all_gene_list.txt')

#Keep only strains that are counted in all datasets.
adip = adip[adip$Strain %in% liver$Strain,]
liver = liver[liver$Strain %in% adip$Strain,]
musc = musc[musc$Strain %in% liver$Strain,]
int = int[int$Strain %in% liver$Strain,]
hypo = hypo[hypo$Strain %in% liver$Strain,]

#Add tissue to gene annotation to later merge data.
adip$gene = paste0(adip$gene_symbol, '_adipose')
liver$gene = paste0(liver$gene_symbol, '_liver')
musc$gene = paste0(musc$Symbol, '_muscle')
int$gene = paste0(int$Symbol, '_intestine')
hypo$gene = paste0(hypo$gene_short_name, '_hypothalamus')

#Simplify the dataframes.
a1 = adip %>% select(Strain, gene, gene_symbol, expression_value)
l1 = liver %>% select(Strain, gene, gene_symbol, expression_value)
m1 = musc %>% select(Strain, gene, 'gene_symbol' = Symbol, expression_value)
i1 = int %>% select(Strain, gene, 'gene_symbol' = Symbol, expression_value)
hypo1 = hypo %>% select(Strain, gene, 'gene_symbol' = gene_short_name, 'expression_value' = value)


#Begin comparing internal correlation structure of gene list - here we will look at all non-circadian genes, listed in #h_pathway$all.liver.NON.circadian.genes
liv_pathway = l1[l1$gene_symbol %in% h_pathway$all.liver.NON.circadian.genes,]
h2 = dcast(liv_pathway, Strain ~ gene_symbol, fun.aggregate = mean, value.var = 'expression_value')
row.names(h2) = h2$Strain
h2$Strain = NULL

#Generate new matrix of midweight bicorrelation coefficients and corresponding statistics
cc1 = bicorAndPvalue(h2, h2, use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p
tt3 = ifelse(tt4 < 0.01,"*","")

#Visualize heatmap of internal correlation structure of all non-circadian genes in liver
pheatmap(cc2, display_numbers = tt3, number_color = "black", color = viridis(50), main='Internal Correlation Structure Non-Circadian Genes', fontsize_row = 3, fontsize_col = 3)

#Print mean cor coef and pvalue
mean(cc2)
mean(tt4)

#Bind all peripheral tissues together and generate dataframe where strains match the non-circadian pathways
all_tissues = as.data.frame(rbind(a1, m1, i1, hypo1))
bin_mat = dcast(all_tissues, Strain ~ gene, value.var = 'expression_value', fun.aggregate = mean)
row.names(bin_mat) = bin_mat$Strain
bin_mat$Strain = NULL
bin_mat = bin_mat[row.names(bin_mat) %in% row.names(h2),]

#Generate the same bicorrelation and pvalue matrices for all peripheral genes : liver non-circadian genes
#This will take a fair amount of time depending on CPU
target_cor = bicorAndPvalue(bin_mat, h2, use = 'p') 

#Extract pvalue matric and calculate -log10(rowMeans), refered to as Ssec score
tt1 = target_cor$p
tt1[is.na(tt1)] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))
colnames(cc3) = 'Ssec_score'
cc3$gene_symbol = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)
scores=cc3
scores$tissue =gsub(".*_","",scores$gene_symbol)
scores$gene_only =gsub("_.*","",scores$gene_symbol)

#Normalize each Ssec score to the mean correlation within each tissue
ss1 = scores
ss2 = ss1[grepl('adipose', ss1$tissue),]
ss2$fitted_ssec = ss2$Ssec_score/mean(ss2$Ssec_score)
ss4 = ss1[grepl('hypothalamus', ss1$tissue),]
ss4$fitted_ssec = ss4$Ssec_score/mean(ss4$Ssec_score)
ss5 = ss1[grepl('intestine', ss1$tissue),]
ss5$fitted_ssec = ss5$Ssec_score/mean(ss5$Ssec_score)
ss6 = ss1[grepl('muscle', ss1$tissue),]
ss6$fitted_ssec = ss6$Ssec_score/mean(ss6$Ssec_score)
new_scores = as.data.frame(rbind(ss2, ss4, ss5, ss6))
merged_scores = new_scores

#Next the analysis repeated for the liver non-restored genes, in order to generate a lists of cross-tissue correlations enriched specifically for this set
#h_pathway$not.restored
liv_pathway = l1[l1$gene_symbol %in% h_pathway$not.restored,]
h2 = dcast(liv_pathway, Strain ~ gene_symbol, fun.aggregate = mean, value.var = 'expression_value')
row.names(h2) = h2$Strain
h2$Strain = NULL

cc1 = bicorAndPvalue(h2, h2, use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p
tt3 = ifelse(tt4 < 0.01,"*","")

pheatmap(cc2, display_numbers = tt3, number_color = "black", color = viridis(50), main='Internal Correlation Structure Not-restored Circadian Genes', fontsize_row = 3, fontsize_col = 3)

mean(cc2)
mean(tt4)

bin_mat = bin_mat[row.names(bin_mat) %in% row.names(h2),]


target_cor = bicorAndPvalue(bin_mat, h2, use = 'p')
tt1 = target_cor$p
tt1[is.na(tt1)] = 0.5
cc3 = as.data.frame(rowMeans(-log10(tt1)))
colnames(cc3) = 'Ssec_score'
cc3$gene_symbol = row.names(cc3)
cc3 = cc3[order(cc3$Ssec_score, decreasing = T),]
summary(cc3$Ssec_score)
scores=cc3
mean(scores$Ssec_score + (sd(scores$Ssec_score)*2))
scores$tissue =gsub(".*_","",scores$gene_symbol)
scores$gene_only =gsub("_.*","",scores$gene_symbol)

ss1 = scores
table(ss1$tissue)
ss2 = ss1[grepl('adipose', ss1$tissue),]
ss2$fitted_ssec = ss2$Ssec_score/mean(ss2$Ssec_score)
ss4 = ss1[grepl('hypothalamus', ss1$tissue),]
ss4$fitted_ssec = ss4$Ssec_score/mean(ss4$Ssec_score)
ss5 = ss1[grepl('intestine', ss1$tissue),]
ss5$fitted_ssec = ss5$Ssec_score/mean(ss5$Ssec_score)
ss6 = ss1[grepl('muscle', ss1$tissue),]
ss6$fitted_ssec = ss6$Ssec_score/mean(ss6$Ssec_score)
new_scores = as.data.frame(rbind(ss2, ss4, ss5, ss6))

#Now the nonrestored genes will be added to the final list and we will calculate fold-enrichment for each peripheral gene (nonrestored Ssec/ non-circadian Ssec)
merged_scores$non_restored_ssec = new_scores$fitted_ssec[match(merged_scores$gene_symbol, new_scores$gene_symbol)]
merged_scores$normalized_nonrestored = merged_scores$non_restored_ssec / merged_scores$fitted_ssec
merged_scores$secreted = match(merged_scores$gene_only, sec_prots$Gene.names...primary.., nomatch=0)
merged_scores$secreted = merged_scores$secreted > 0
scores = merged_scores
scores = scores[order(scores$normalized_nonrestored, decreasing=T),]
ssec_table = scores

#Write this table as a list
write.table(ssec_table, file = 'Cross-tissue correlations with nonrestored liver clock genes - normalized to global non-circadian ssecs.txt', sep = '\t', row.names=F)

#Generate color scheme for plotting
colors <- c("gold", "red", "dodgerblue1","green")
names(colors) = c("adipose", "muscle", "intestine", "hypothalamus")
scores$tissue_col = colors[match(scores$tissue, names(colors))]

#Inspect the top genes
head(scores)

#plot all peripheral gene enrichments for nonrestored liver genes - this plot will take some time to generate
ggplot(scores, aes(x=fct_reorder2(gene_symbol, normalized_nonrestored, normalized_nonrestored, .desc = T), y=normalized_nonrestored)) + geom_col(fill=scores$tissue_col) + theme(axis.text.x = element_text(angle=90, size=8), plot.title = element_text(hjust=0.5)) +ggtitle('All genes enriched for liver not_restored circadian pathways') + theme_minimal() + geom_hline(yintercept=mean(scores$normalized_nonrestored+ (sd(scores$normalized_nonrestored)*2)), linetype="dashed", color = "gray9",cex=3)+ xlab('gene') + ylab('cross-tissue scorec (Ssec) / non-liver ssec') + theme_classic() ### this will take a long time ### 

#take top 20 peripheral genes and plot
scores1 = scores[1:20,]
ggplot(scores1, aes(x=fct_reorder2(gene_symbol, normalized_nonrestored, normalized_nonrestored, .desc = T), y=normalized_nonrestored)) + geom_col(fill=scores1$tissue_col) +  ggtitle('All genes enriched for liver not_restored circadian pathways') + theme_minimal()  + xlab('gene') + ylab('cross-tissue scorec (Ssec) / non-circadian ssec')+ theme(axis.text.x = element_text(angle=90, size=8), plot.title = element_text(hjust=0.5)) + geom_hline(yintercept=mean(scores$normalized_nonrestored + (sd(scores$normalized_nonrestored)*2)), linetype="dashed", color = "gray9",cex=3)

#Set a cutoff (Here we use 3 standard deviations above mean Ssec) and generate pie chart for distribution of highly significant genes
scores1 = scores[scores$normalized_nonrestored > mean(scores$normalized_nonrestored + (sd(scores$normalized_nonrestored)*3)),]
table(scores1$tissue)
pie(table(scores1$tissue))





