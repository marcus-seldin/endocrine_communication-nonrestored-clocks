library(dplyr)
library(reshape2)
library(viridis)
library(WGCNA)
library(pheatmap)
library(ggplot2)
library(forcats)
library(colormap)
allowWGCNAThreads()
sec_prots = read.delim('secreted_proteins.txt')

adip = read.delim('HMDP_HFHS_trx_adipose.txt')
musc = read.delim('HMDP_HFHS_trx_muscle.txt')
liver = read.delim('HMDP_HFHS_trx_liver.txt')
int = read.delim('HMDP_HFHS_trx_intestine.txt')
hypo = read.delim('HMDP_HFHS_trx_hypothalamus.txt')
traits = read.delim('HMDP_HFHS_traits.txt')

adip = adip[adip$Strain %in% liver$Strain,]
liver = liver[liver$Strain %in% adip$Strain,]
musc = musc[musc$Strain %in% liver$Strain,]
int = int[int$Strain %in% liver$Strain,]
hypo = hypo[hypo$Strain %in% liver$Strain,]
trait = traits[traits$Strain %in% liver$Strain,]

adip$gene = paste0(adip$gene_symbol, '_adipose')
liver$gene = paste0(liver$gene_symbol, '_liver')
musc$gene = paste0(musc$Symbol, '_muscle')
int$gene = paste0(int$Symbol, '_intestine')
hypo$gene = paste0(hypo$gene_short_name, '_hypothalamus')

a1 = adip %>% select(Strain, gene, gene_symbol, expression_value)
l1 = liver %>% select(Strain, gene, gene_symbol, expression_value)
m1 = musc %>% select(Strain, gene, 'gene_symbol' = Symbol, expression_value)
#hrt1 = heart %>% select(Strain, gene, gene_symbol, 'expression_value' = kallisto_tpm)
i1 = int %>% select(Strain, gene, 'gene_symbol' = Symbol, expression_value)
hypo1 = hypo %>% select(Strain, gene, 'gene_symbol' = gene_short_name, 'expression_value' = value)


trts = dcast(trait, Strain ~ trait_name, value.var = 'value', fun.aggregate = mean)
h_pathway = read.delim('G:/My Drive/lab files/paolo/autonomous clock lookups/all_gene_list.txt')
l1$gene_symbol[1:10]
h_pathway$all.liver.NON.circadian.genes
heart_pathway = l1[l1$gene_symbol %in% h_pathway$all.liver.NON.circadian.genes,]
head(heart_pathway)
h2 = dcast(heart_pathway, Strain ~ gene_symbol, fun.aggregate = mean, value.var = 'expression_value')
row.names(h2) = h2$Strain
h2$Strain = NULL

cc1 = bicorAndPvalue(h2, h2, use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p
tt3 = ifelse(tt4 < 0.01,"*","")

#pheatmap(cc2, display_numbers = tt3, number_color = "black", color = viridis(50), main='Internal Correlation Structure Non-Circadian Genes', fontsize_row = 3, fontsize_col = 3)

mean(cc2)
mean(tt4)
all_tissues = as.data.frame(rbind(a1, m1, i1, hypo1))
bin_mat = dcast(all_tissues, Strain ~ gene, value.var = 'expression_value', fun.aggregate = mean)
row.names(bin_mat) = bin_mat$Strain
bin_mat$Strain = NULL
row.names(bin_mat)[1:10]
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
#ss3 = ss1[grepl('heart', ss1$tissue),]
#ss3$fitted_ssec = ss3$Ssec_score/mean(ss3$Ssec_score)
ss4 = ss1[grepl('hypothalamus', ss1$tissue),]
ss4$fitted_ssec = ss4$Ssec_score/mean(ss4$Ssec_score)
ss5 = ss1[grepl('intestine', ss1$tissue),]
ss5$fitted_ssec = ss5$Ssec_score/mean(ss5$Ssec_score)
ss6 = ss1[grepl('muscle', ss1$tissue),]
ss6$fitted_ssec = ss6$Ssec_score/mean(ss6$Ssec_score)

new_scores = as.data.frame(rbind(ss2, ss4, ss5, ss6))
merged_scores = new_scores

#h_pathway$not.restored
heart_pathway = l1[l1$gene_symbol %in% h_pathway$not.restored,]
h2 = dcast(heart_pathway, Strain ~ gene_symbol, fun.aggregate = mean, value.var = 'expression_value')
row.names(h2) = h2$Strain
h2$Strain = NULL

cc1 = bicorAndPvalue(h2, h2, use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p
tt3 = ifelse(tt4 < 0.01,"*","")

pheatmap(cc2, display_numbers = tt3, number_color = "black", color = viridis(50), main='Internal Correlation Structure Not-restored Circadian Genes', fontsize_row = 3, fontsize_col = 3)

mean(cc2)
mean(tt4)
#all_tissues = as.data.frame(rbind(a1, m1, i1, hypo1))
#bin_mat = dcast(all_tissues, Strain ~ gene, value.var = 'expression_value', fun.aggregate = mean)
#row.names(bin_mat) = bin_mat$Strain
#bin_mat$Strain = NULL
#row.names(bin_mat)[1:10]

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
merged_scores$non_restored_ssec = new_scores$fitted_ssec[match(merged_scores$gene_symbol, new_scores$gene_symbol)]
merged_scores$normalized_nonrestored = merged_scores$non_restored_ssec / merged_scores$fitted_ssec
merged_scores$secreted = match(merged_scores$gene_only, sec_prots$Gene.names...primary.., nomatch=0)
merged_scores$secreted = merged_scores$secreted > 0


scores = merged_scores
scores = scores[order(scores$normalized_nonrestored, decreasing=T),]
ssec_table = scores
write.table(ssec_table, file = 'Cross-tissue correlations with nonrestored liver clock genes - normalized to global non-circadian ssecs.txt', sep = '\t', row.names=F)

colors <- c("darkorange1", "darkorchid4","dodgerblue1","deeppink1")
names(colors) = c("adipose", "muscle", "intestine", "hypothalamus")
scores$tissue_col = colors[match(scores$tissue, names(colors))]

head(scores)
#ggplot(scores, aes(x=fct_reorder2(gene_symbol, normalized_nonrestored, normalized_nonrestored, .desc = T), y=normalized_nonrestored)) + geom_col(fill=scores$tissue_col) + theme(axis.text.x = element_text(angle=90, size=8), plot.title = element_text(hjust=0.5)) +ggtitle('All genes enriched for liver not_restored circadian pathways') + theme_minimal() + geom_hline(yintercept=mean(scores$normalized_nonrestored+ (sd(scores$normalized_nonrestored)*2)), linetype="dashed", color = "gray9",cex=3)+ xlab('gene') + ylab('cross-tissue scorec (Ssec) / non-liver ssec') + theme_classic()


scores1 = scores[1:20,]
ggplot(scores1, aes(x=fct_reorder2(gene_symbol, normalized_nonrestored, normalized_nonrestored, .desc = T), y=normalized_nonrestored)) + geom_col(fill=scores1$tissue_col) +  ggtitle('All genes enriched for liver not_restored circadian pathways') + theme_minimal()  + xlab('gene') + ylab('cross-tissue scorec (Ssec) / non-circadian ssec')+ theme(axis.text.x = element_text(angle=90, size=8), plot.title = element_text(hjust=0.5)) + geom_hline(yintercept=mean(scores$normalized_nonrestored + (sd(scores$normalized_nonrestored)*2)), linetype="dashed", color = "gray9",cex=3)


#pie chart for distribution of highly significant genes
scores1 = scores[scores$normalized_nonrestored > mean(scores$normalized_nonrestored + (sd(scores$normalized_nonrestored)*3)),]
table(scores1$tissue)
pie(table(scores1$tissue))
#Look at circadian vs non-circadian scores in not-restored
#scores1 = scores[scores]


bin_mat = dcast(all_tissues, Strain ~ gene, value.var = 'expression_value', fun.aggregate = mean)
row.names(bin_mat) = bin_mat$Strain
bin_mat$Strain = NULL
row.names(bin_mat)[1:10]
bin_mat = bin_mat[row.names(bin_mat) %in% row.names(h2),]
pp1 = bicorAndPvalue(bin_mat, h2, use='p')
nonr_cors = melt(pp1$p)
head(nonr_cors)
colnames(nonr_cors) = c('origin_gene', 'liver_gene', 'pvalue')
nonr_cors$log10p = -log10(nonr_cors$pvalue)
nonr_cors$origin_gene_only = gsub("_.*","", nonr_cors$origin_gene) 
nonr_cors$pathway = match(nonr_cors$origin_gene_only, h_pathway$all.liver.circadian.genes, nomatch = 0)
nonr_cors$pathway = ifelse(nonr_cors$pathway > 0, 'circadian', 'non_circadian')
nonr_cors = na.omit(nonr_cors)
tt1 = nonr_cors %>% select(pathway, log10p) %>% group_by(pathway) %>% summarize(avg=mean(log10p), n=n(), sd = sd(log10p), se=sd/sqrt(n), lower.ci = avg - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = avg + qt(1 - (0.05 / 2), n - 1) * se)
tt1$pathway = factor(tt1$pathway, levels=c('non_circadian', 'circadian'))
#maybe count the numebr of significant??
ggplot(tt1, aes(x=pathway, y=avg, fill=pathway)) + 
  geom_point(stat="identity", color='black', size=3) +
  geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.4,
                position=position_dodge(.9)) + theme_minimal() + ylab('log10(gene-gene bicor pvalue)') +xlab('none')

ggplot(nonr_cors, aes(x=pathway, y=log10p, fill=pathway))  + geom_violin() + geom_boxplot(aes(fill = pathway), width=0.1)






