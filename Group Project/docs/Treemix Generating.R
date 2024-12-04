library(vcf2others)
library(ape)
metatree <- metadata[metadata$id %in% colnames(final_merged.vcf@gt[,-1]),]

vcf2treemix(
  final_merged.vcf,
  ind_pop = ind_pop,
  keep_pop = keepers,
  inc_missing = TRUE,
  out_file = "~/projects/eco_genomics/Group Project/outputs/treemixMerged_infile.txt"
)

treemixMerged <- read.table("~/projects/eco_genomics/Group Project/outputs/treemixMerged_infile.txt")

Phylo_data <- read.tree(text = "(((SEU:0.00504761,CEU:0.00151707):0.00153618,PNW:0.00238455):0.000119893,(OG:0.101573,((NE:0.000413659,NEU:0.00306033):0.001276,WEU:0.000759046):0.000745533):0.000119
893)")
plot_tree("out_stem")
plot_tree("out_stem3")
plot_tree("out_stem4")
#No migration events tree: ln(likelihood) = 167.83304
#One migration events tree: ln(likelihood) = 194.99464
#Two migration events tree: ln(likelihood) = 199.97198
#Three migration events tree: ln(likelihood) = 200.92086
#Four migration events tree: ln(likelihood) = 201.48137
vcf2treemix(
  final_merged.vcf,
  ind_pop = ind_pop,
  keep_pop = keepers,
  inc_missing = FALSE,
  out_file = "~/projects/eco_genomics/Group Project/outputs/treemixMergedNoMissing_infile.txt"
)
TreemixMergedNoMissing <- read.table("treemixMergedNoMissing_infile.txt")
