

#---Load packages
library(here)
library(ggplot2)
library(tidyverse)
library(ggpubr)

# Define custom palette
cpal <- c("darkgoldenrod3", "dimgrey", "darkgreen")

#----Plot number of SNPs per species----

# Load atlas metadata info
load(here("data", "stress_samples_for_DGE_analysis.RData"))
stress_info <- Reduce(rbind, biotic_samplelist)
stress_info$Pathogen[stress_info$Pathogen == "SBA"] <- "Aglycines"
stress_info$Pathogen[stress_info$Pathogen == "Pgregata"] <- "Cgregata"

# Load SNP information
load(here("products", "result_files", "snp_granges.rda"))
snp_info <- as.data.frame(snp_granges)


# Create a stress class ontology for each pest
stress_ontology <- data.frame(
    Pathogen = c("Aglycines", "SBA", "Fgraminearum", "Foxysporum",
                 "Fvirguliforme", "Hglycines", "MAMP", "Mphaseolina",
                 "Pgregata", "Ppachyrhizi", "Psojae", "Psojae_glucan_elicitor",
                 "Rreniformis", "Rsolani", "Slitura", "SMV",
                 "Ssclerotiorum", "Evarivestis", "Efabae", "Pincludens",
                 "Agemmatalis", "Xaxonopodis", "Cgregata", "Dphaseolorum",
                 "BPMV", "PMV", "Mincognita", "TRSV", 
                 "Psylvaticum", "Fequiseti"),
    Class = c("insect", "insect", "fungus", "fungus",
              "fungus", "nematode", "MAMP", "fungus", 
              "fungus", "fungus", "oomycete", "oomycete",
              "nematode", "fungus", "insect", "virus",
              "fungus", "insect", "insect", "insect",
              "insect", "bacterium", "fungus", "fungus",
              "virus", "virus", "nematode", "virus",
              "oomycete", "fungus")
)

# Count occurrences for SNP and transcriptome
st1 <- stress_ontology %>%
    inner_join(snp_info, by = c("Pathogen" = "Organism")) %>%
    filter(Class == "insect" | Class == "nematode") %>%
    distinct(SNP, .keep_all=TRUE)
write_tsv(st1, here("products", "tables", "sup_table1.tsv"))

count_snp <- stress_ontology %>%
    inner_join(snp_info, by = c("Pathogen" = "Organism")) %>%
    distinct(SNP, .keep_all=TRUE) %>%
    filter(Class == "insect" | Class == "nematode") %>%
    dplyr::select(Pathogen, Class) %>%
    group_by(Pathogen) %>%
    summarise(n = n())

st2 <- stress_ontology %>%
    inner_join(stress_info)
write_tsv(st2, here("products", "tables", "sup_table2.tsv"))

count_transcriptome <- stress_ontology %>%
    inner_join(stress_info) %>%
    filter(Class == "insect" | Class == "nematode") %>%
    dplyr::select(Pathogen, Class) %>%
    replace("SBA", "Aglycines") %>%
    group_by(Pathogen) %>%
    summarise(n = n())

# Merge data frames by pathogen 
count_snp_transcriptome <- full_join(count_snp, count_transcriptome,
                                     by = c("Pathogen" = "Pathogen")) %>%
    dplyr::rename(SNP = n.x, 
           Transcriptome = n.y) %>%
    replace(is.na(.), 0) %>%
    filter(SNP > 5 | Transcriptome > 5) %>%
    left_join(stress_ontology)

# Plot data for pests
pest <- count_snp_transcriptome %>%
    pivot_longer(!c(Pathogen, Class)) %>%
    filter(Class == "insect" | Class == "nematode") %>%
    arrange(Pathogen) %>%
    mutate(Pathogen = fct_relevel(Pathogen)) %>%
    dplyr::rename(Species = Pathogen, Frequency = value) %>%
    ggbarplot(., x="Species", y="Frequency", fill="Species",
              palette=ggsci::pal_d3("category20")(20),
              facet.by="name", orientation="horiz",
              legend="none")

pest_overlap <- count_snp_transcriptome %>%
    filter(SNP > 0 & Transcriptome > 0) %>%
    pivot_longer(!c(Pathogen, Class)) %>%
    filter(Class == "insect" | Class == "nematode") %>%
    arrange(Pathogen) %>%
    mutate(Pathogen = fct_relevel(Pathogen)) %>%
    dplyr::rename(Species = Pathogen, Frequency = value) %>%
    ggbarplot(., x="Species", y="Frequency", fill="Species", color="Species",
              palette = cpal,
              facet.by="name", orientation="horiz",
              legend="none") +
    theme_bw() +
    ggtitle("SNPs and RNA-seq samples") +
    theme(plot.title = element_text(face="bold", size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")

ggsave(insect_overlap, width=6, height=6,
       filename = here("products", "plots", 
                       "frequency_of_snps_and_transcriptome_samples_overlap.png"))

#----Plot SNP positions in the genome----
library(GenomicRanges)
library(GenomeInfoDb)

# Load SNP ranges
load(here("products", "result_files", "snp_granges.rda"))

# Load gff and remove scaffolds
R.utils::gunzip(here("data", "PLAZA_selected.transcripts.gff.gz"), remove=FALSE)
gff <- gread::read_gff(here("data", "PLAZA_selected.transcripts.gff"))
file.remove(here("data", "PLAZA_selected.transcripts.gff"))
scaffolds <- grep("scaffold", seqlevels(gff), value=TRUE)
gff <- dropSeqlevels(gff, scaffolds, pruning.mode = "tidy")

# Create a different object for each genome location
table(gff$feature)


granges_all <- gread::construct_introns(gff, update=TRUE)
introns <- gread::construct_introns(gff, update=FALSE)
exons <- gff[gff$feature == "exon", ]
fivep_utr <- gff[gff$feature == "five_prime_UTR", ]
threep_utr <- gff[gff$feature == "three_prime_UTR", ]

snp_grangeslist <- lapply(snp_grangeslist, function(x) {
    y <- x
    end(y) <- start(y)
    return(y)
})

# Create a data frame containing the position of each SNP in the genome
# Intron
snps_intron <- sapply(snp_grangeslist, function(x) {
    y <- IRanges::subsetByOverlaps(x, introns)$SNP
    return(y)
})
n_intron <- sapply(snps_intron, length)
snps_intron_df <- data.frame(
    SNPs = unlist(snps_intron),
    Location = "intron",
    Species = c(rep("Aglycines", n_intron[1]),
                rep("Hglycines", n_intron[2]),
                rep("Slitura", n_intron[3]))
)

# Exon
snps_exon <- sapply(snp_grangeslist, function(x) {
    y <- IRanges::subsetByOverlaps(x, exons)$SNP
    return(y)
})
n_exon <- sapply(snps_exon, length)
snps_exon_df <- data.frame(
    SNPs = unlist(snps_exon),
    Location = "exon",
    Species = c(rep("Aglycines", n_exon[1]),
                rep("Hglycines", n_exon[2]),
                rep("Slitura", n_exon[3]))
)

# 5'-UTR
snps_5p <- sapply(snp_grangeslist, function(x) {
    y <- IRanges::subsetByOverlaps(x, fivep_utr)$SNP
    return(y)
})
n_5p <- sapply(snps_5p, length)
snps_5p_df <- data.frame(
    SNPs = unlist(snps_5p),
    Location = "5p",
    Species = c(rep("Aglycines", n_5p[1]),
                rep("Hglycines", n_5p[2]),
                rep("Slitura", n_5p[3]))
)

# 3'-UTR
snps_3p <- sapply(snp_grangeslist, function(x) {
    y <- IRanges::subsetByOverlaps(x, threep_utr)$SNP
    return(y)
})
n_3p <- sapply(snps_3p, length)
snps_3p_df <- data.frame(
    SNPs = unlist(snps_3p),
    Location = "3p",
    Species = c(rep("Aglycines", n_3p[1]),
                rep("Hglycines", n_3p[2]),
                rep("Slitura", n_3p[3]))
)

# Intergenic
snps_intergenic <- sapply(seq_along(snp_grangeslist), function(x) {
    others <- c(snps_intron[[x]], snps_exon[[x]], snps_5p[[x]], snps_3p[[x]])
    y <- snp_grangeslist[[x]]$SNP[!snp_grangeslist[[x]]$SNP %in% others]
    return(y)
})
n_intergenic <- sapply(snps_intergenic, length)
snps_intergenic_df <- data.frame(
    SNPs = unlist(snps_intergenic),
    Location = "Intergenic",
    Species = c(rep("Aglycines", n_intergenic[1]),
                rep("Hglycines", n_intergenic[2]),
                rep("Slitura", n_intergenic[3]))
)

snps_and_location <- rbind(
    snps_intron_df, snps_exon_df, snps_5p_df, snps_3p_df, snps_intergenic_df
)
readr::write_tsv(
    snps_and_location, 
    file = here("products", "result_files", "snps_and_location.tsv")
)


location_of_snps <- data.frame(
    Intergenic = sapply(snps_intergenic, length),
    Intron = sapply(snps_intron, length),
    Exon = sapply(snps_exon, length),
    Fiveprime_UTR = sapply(snps_5p, length),
    Threeprime_UTR = sapply(snps_3p, length),
    Organism = names(snp_grangeslist)
)

head(location_of_snps)

# Melt data frame
location_of_snps_plotdata <- reshape2::melt(location_of_snps)
location_of_snps_plotdata$variable <- factor(location_of_snps_plotdata$variable, 
                                             levels=c("Intergenic", "Exon", 
                                                      "Intron",
                                                      "Fiveprime_UTR", 
                                                      "Threeprime_UTR"))

snp_positions <- location_of_snps_plotdata %>%
    ggpubr::ggbarplot(data=.,
                      x="variable", y="value", 
                      facet.by = "Organism",
                      fill="Organism", color="Organism", 
                      palette = cpal,
                      orientation="horiz", ncol=5, legend="none",
                      xlab="Genome location", ylab="Absolute frequency",
                      label=TRUE, lab.pos="out", lab.hjust = -0.3, lab.vjust = 0.1,
                      title="Location of SNPs in the genome", font.title="bold") +
    ggplot2::expand_limits(y = 70) + 
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::theme(axis.text.x = element_text(size=11)) +
    scale_x_discrete(labels=c("Threeprime_UTR" = "3'-UTR", 
                              "Fiveprime_UTR" = "5'-UTR",
                              "Intron" = "Intron",
                              "Exon" = "Exon",
                              "Intergenic" = "Intergenic")) +
    theme_bw() +
    theme(plot.title = element_text(face="bold", size=12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(angle=45, hjust=1))
    

snp_positions


#----Circos plot of SNP positions across chromosomes per stress class----
library(cageminer)
load(here("data", "soybean_genome_ranges.rda"))
gff_gene <- gff[gff$feature == "gene", ]
snp_grangeslist <- GRangesList(snp_grangeslist)
snp_circos <- plot_snp_circos(genome.ranges, gff_gene, snp_grangeslist) +
    ggplot2::scale_color_manual(values = cpal)

#----Plot SNP distribution across chromosomes----
library(ggsci)
snp_dist <- plot_snp_distribution(snp_grangeslist) +
    scale_x_discrete(labels = rev(c(
        "Gm20", "Gm19", "Gm18", "Gm17", "Gm16", 
        "Gm15", "Gm14", "Gm13", "Gm12", "Gm11",
        "Gm10", "Gm09", "Gm08", "Gm07", "Gm06", 
        "Gm05", "Gm04", "Gm03", "Gm02", "Gm01" 
    ))) +
    scale_fill_manual(values = cpal)

# Save plots individually
snp_circos # Save manually in RStudio with width=800 and height=800
ggsave(snp_dist, filename = here("products", "plots", "snp_dist.png"),
       width=8, height=4)
ggsave(snp_positions, filename = here("products", "plots", "snp_positions.png"),
       width=8, height=2.5)


library(ggpubr)
library(magick)
circos <- magick::image_read(here("products", "plots", "snp_circos.png"))
circosgg <- magick::image_ggplot(circos, interpolate = TRUE)

panel1 <- ggarrange(pest_overlap, circosgg, ncol=2, widths = c(2,3),
                    labels = c("A", "B"))
panel2 <- ggarrange(snp_dist, snp_positions, nrow=2, heights = c(1.5, 1),
                    labels=c("C", "D"))
final_plot <- ggarrange(panel1, panel2, heights = c(1, 1.5), nrow=2)
ggsave(final_plot, 
       filename = here("products", "plots", "main_stats_SNPs.pdf"),
       width=9, height=12)










