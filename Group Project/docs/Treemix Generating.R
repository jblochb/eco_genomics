
############ This script was used only for the treemix analyses associated with the merged vcf that would provide an out group on the tree. It includes the steps to merge the vcfs. 
#We did not end up using this script in our final paper because we are not confident in the migration events being inferred between the out group and the regions. 


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

treemixMerged <- read.table("~/projects/eco_genomics/Group Project/outputs/treemixMerged_infile.txt.gz")

Phylo_data <- read.tree(text = "(((SEU:0.00504761,CEU:0.00151707):0.00153618,PNW:0.00238455):0.000119893,(OG:0.101573,((NE:0.000413659,NEU:0.00306033):0.001276,WEU:0.000759046):0.000745533):0.000119
893)")
plot_tree("out_stem")
plot_tree("out_stem3")
plot_tree("out_stem4")
plot_tree("out_stem5")
#No migration events tree: ln(likelihood) = 167.83304
#One migration events tree: ln(likelihood) = 194.99464
#Two migration events tree: ln(likelihood) = 199.97198
#Three migration events tree: ln(likelihood) = 200.92086
#Four migration events tree: ln(likelihood) = 201.48137
#Five migration events tree: ln(likelihood) = 202.15676
vcf2treemix(
  final_merged.vcf,
  ind_pop = ind_pop,
  keep_pop = keepers,
  inc_missing = FALSE,
  out_file = "~/projects/eco_genomics/Group Project/outputs/treemixMergedNoMissing_infile.txt"
)
TreemixMergedNoMissing <- read.table("treemixMergedNoMissing_infile.txt")


#################  Filtering the ref vcf for better outgroup data

ref_two_allele_filt.vcf <- filter_biallelic(ref.vcf)
#min_mac = minimum minor allele count (this gets rid of NAs)
ref_two_allele_filt.vcf <- min_mac(ref_two_allele_filt.vcf, min.mac = 1)
ref_two_allele_filt.vcf
ref.SNPMiss.vcf <- missing_by_snp(ref_two_allele_filt.vcf, cutoff = 0.5)


################################## Extracting common SNPs

gt_vcf.thin <- extract.gt(vcf.thin)
gt_ref.vcf.thin <- extract.gt(reference.vcf)

# Find common SNPs between the two VCF files
common_snps <- intersect(rownames(gt_vcf.thin), rownames(gt_ref.vcf.thin))

# Subset the genotype data to include only the common SNPs
gt1_common <- gt_vcf.thin[common_snps, ]
gt2_common <- gt_ref.vcf.thin[common_snps, ]

# Combine the common SNP genotypes
merged_gt_common <- cbind(gt1_common, gt2_common)

#Create snp ids csv
snp_id_df <- merged_gt_common[,1]
snp_id_df <- as.data.frame(snp_id_df)
snp_ids <- as.character(rownames(snp_id_df))
snp_ids <- as.data.frame(snp_ids)
view(snp_ids)
write.table(snp_ids, "~/projects/eco_genomics/Group Project/outputs/common.pos.txt", quote = F, row.names = F)
write_csv(snp_ids, "~/projects/eco_genomics/Group Project/outputs/snp_ids")
write.vcf(ref_two_allele_filt.vcf, "~/projects/eco_genomics/Group Project/outputs/reference.vcf.gz")
reference.vcf <- read.vcfR("~/projects/eco_genomics/Group Project/outputs/reference.vcf.gz")
merged_outgroup.vcf <- read.vcfR("~/projects/eco_genomics/Group Project/outputs/merged_outgroup.vcf.gz.recode.vcf")

metatree <- metadata[metadata$id %in% colnames(merged_outgroup.vcf@gt[,-1]),]

vcf2treemix(
  merged_outgroup.vcf,
  ind_pop = ind_pop,
  keep_pop = keepers,
  inc_missing = TRUE,
  out_file = "~/projects/eco_genomics/Group Project/outputs/treemix170_infile.txt"
)
