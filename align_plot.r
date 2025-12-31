

setwd("/share/nas5/huangls/test/genome_denovo/S288c/40.assembly_stat/")

library(pafr, quietly=TRUE)

ali <- read_paf("next.paf")

#ali <- subset(ali, alen > 1e4)



dotplot(ali)

mean(ali$alen)

ggplot(ali, aes(alen)) + 
  geom_histogram(colour="black", fill="steelblue", bins=20) + 
  theme_bw(base_size=16) + 
  ggtitle("Distribution of alignment lengths") +
  scale_x_log10("Alignment-length")


long_ali <- subset(ali, alen > 1e4)
plot_synteny(long_ali, q_chrom=ali$qname, t_chrom=ali$tname, centre=TRUE)


plot_coverage(long_ali, fill='qname') 
  #scale_fill_brewer(palette="Set1")