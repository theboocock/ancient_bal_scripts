library(foreach)
library(SNPRelate)
library(ape)
library(stringr)
library(phangorn)
library(ggtree)
library(treeio)
library(ggplot2)
library(tidyverse)
library(glue)
genofile <- snpgdsOpen("data/vcf//filtered.gds")
chinese_annotations1 = read.csv("data/annotations/chinese_annotations.csv",header=T)
chinese_annotations2 = read.delim("data/annotations/chinese_annotations_2.txt", header=T, sep="\t")
chinese_annotations_merged = merge(chinese_annotations1,chinese_annotations2, by.x="BioSample.ID",by.y="BioSample")
chinese_annotations_merged$name = substr(chinese_annotations_merged$Genome.Accession,0,4)
annotation_gal_genes = read.csv("data/annotations/gal_annotations_cbs_complete.csv",stringsAsFactors = F)
annotation_gal_genes$GAL2[annotation_gal_genes$gal2 == "YLR081W"] = "wt"
colnames(annotation_gal_genes) = c("Standardized.name","GAL1","GAL10","GAL7","GAL2","PGM1")
annotation_gal_genes$alt_id = annotation_gal_genes$Standardized.name
annotation_gal_genes$alt_id[na.omit(match(chinese_annotations_merged$name,annotation_gal_genes$Standardized.name))] = as.character(chinese_annotations_merged$Run)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
calculate.pi = function(gt){
  n = nrow(gt)
  average_distance = sum(ape::dist.gene(gt))/(n*(n - 1)/2)
  return(average_distance)
  
}
calculate.random_pi = function(vcf_in, strains, number=100, length=7296,genome="data/vcf/chrom_lengths.txt"){
  pis = c()
  command_str = glue("bedtools random -n {number} -l {length}  -g {genome}",number=number, length=length, genome=genome)
  tabix_regions = system(command_str,intern = T)
  output_df = data.frame()
  for(tabix_r in tabix_regions){
    split_tmp = str_split(tabix_r,"\t")[[1]]
    chrom = split_tmp[1]
    start = split_tmp[2]
    end =  split_tmp[3]
    tabix_string_for_command = glue("{chrom}:{start}-{end}",chrom=chrom,start=start,end=end)
    tabix_command = glue("tabix -h {vcf} {region}", vcf=vcf.fn, region=tabix_string_for_command)
    system(glue("{tabix_command} > tmp/tmp.vcf"))
    snpgdsVCF2GDS("tmp/tmp.vcf", "tmp/tmp5.gds", method="biallelic.only")
    genofile.2 <- snpgdsOpen("tmp/tmp5.gds")
    if(length(read.gdsn(index.gdsn(genofile.2,"snp.position"))) != 0){
      gt = read.gdsn(index.gdsn(genofile.2,"genotype"))
      gt_gal = gt[which(strain_names_vcf %in% strains),]
      pi = calculate.pi(gt_gal)
      print(pi)
      pis = c(pi, pis)
      df_tmp = data.frame(chrom=chrom,start=start,end=end, pi=pi)
      output_df = rbind(output_df,df_tmp)
      
    }
    snpgdsClose(genofile.2)
  }
  return(output_df)
}

gal_strains = annotation_gal_genes$alt_id[annotation_gal_genes$GAL1 == "cbs"]

#tabix_string_for_command = glue("{chrom}:{start}-{end}",chrom=chrom,start=start,end=end)
tabix_command = glue("tabix -h {vcf} alt_gal10_7", vcf=vcf.fn, region=tabix_string_for_command)
system(glue("{tabix_command} > tmp/tmp.vcf"))
snpgdsVCF2GDS("tmp/tmp.vcf", "tmp/tmp.gds", method="biallelic.only")
genofile.2 <- snpgdsOpen("tmp/tmp.gds")
strain_names = read.gdsn(index.gdsn(genofile.2,"sample.id"))
strain_names_vcf = str_replace(strain_names,pattern = "SACE_",replacement = "")
gt = read.gdsn(index.gdsn(genofile.2,"genotype"))
gt_gal = gt[which(strain_names_vcf %in% gal_strains),]
pi_gal = calculate.pi(gt_gal)
snpgdsClose(genofile.2)
#pis_df$pi2 = pis_df$pis / 7296
#pi_gal2 = pi_gal/7296
pis_gal1_10_7 = calculate.random_pi(vcf_in = vcf.fn,strains=gal_strains,length=7296)

#pis_df = data.frame(pis=pis)
#pis_df %>% ggplot(aes(x=pi2)) + geom_histogram() + xlim(c(0,200/7296)) + geom_vline(xintercept = pi_gal/7296, col="red") + theme_bw() + xlab("Nucleotide diversity") +ylab("Count") +
#  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))

### GAL7 promoter

pgm1_strains = annotation_gal_genes$alt_id[annotation_gal_genes$PGM1 == "cbs"]
#abix_string_for_command = glue("{chrom}:{start}-{end}",chrom=chrom,start=start,end=end)
tabix_command = glue("tabix -h {vcf} alt_gal10_7", vcf=vcf.fn, region=tabix_string_for_command)
system(glue("{tabix_command} > tmp/tmp.vcf"))
snpgdsVCF2GDS("tmp/tmp.vcf", "tmp/tmp3.gds", method="biallelic.only")
genofile.2 <- snpgdsOpen("tmp/tmp3.gds")
gt = read.gdsn(index.gdsn(genofile.2,"genotype"))
gal7_1_10_length = 7296
gt_gal = gt[which(strain_names_vcf %in% pgm1_strains),]
pi_mini_gal7 = calculate.pi(gt_gal)
snpgdsClose(genofile.2)
pis_mini_gal7 = calculate.random_pi(vcf_in = vcf.fn,strains=pgm1_strains,length=gal7_1_10_length)

#write.table("")
gal2_strains = annotation_gal_genes$alt_id[(annotation_gal_genes$GAL2 == "cbs")]
#tabix_string_for_command = glue("{chrom}:{start}-{end}",chrom=chrom,start=start,end=end)
tabix_command = glue("tabix -h {vcf} alt_gal2", vcf=vcf.fn, region=tabix_string_for_command)
system(glue("{tabix_command} > tmp/tmp.vcf"))
#system(('cat tmp.vcf | sed "s/\t1|/\t1\//g" | sed "s/\t0|/\t0\//g" > tmp2.vcf'))
snpgdsVCF2GDS("tmp/tmp.vcf", "tmp/tmp.gds", method="biallelic.only")

#snpgdsVCF2GDS("tmp/ars.vcf", "tmp/tmp2.gds", method="biallelic.only")
genofile.2 <- snpgdsOpen("tmp/tmp.gds")
gt = read.gdsn(index.gdsn(genofile.2,"genotype"))

gal2_length = 5682
gt_gal = gt[which(strain_names_vcf %in% gal2_strains),]
pi_gal2 = calculate.pi(gt_gal)
snpgdsClose(genofile.2)

pis_gal2 = calculate.random_pi(vcf_in = vcf.fn,strains=gal2_strains,length=gal2_length)

#pis_pgm1 = calculate.random_pi(vcf_in = vcf.fn,strains=pgm1_strains,length=pgm1_leng
pgm1_strains = annotation_gal_genes$alt_id[annotation_gal_genes$PGM1 == "cbs"]
#abix_string_for_command = glue("{chrom}:{start}-{end}",chrom=chrom,start=start,end=end)
tabix_command = glue("tabix -h {vcf} alt_pgm1", vcf=vcf.fn, region=tabix_string_for_command)
system(glue("{tabix_command} > tmp/tmp.vcf"))
snpgdsVCF2GDS("tmp/tmp.vcf", "tmp/tmp3.gds", method="biallelic.only")
genofile.2 <- snpgdsOpen("tmp/tmp3.gds")
gt = read.gdsn(index.gdsn(genofile.2,"genotype"))
pgm1_length = 2956
gt_gal = gt[which(strain_names_vcf %in% pgm1_strains),]
pi_pgm1 = calculate.pi(gt_gal)
snpgdsClose(genofile.2)
pis_pgm1 = calculate.random_pi(vcf_in = vcf.fn,strains=pgm1_strains,length=pgm1_length,number = 10)
pis_df = data.frame(gal2=pis_gal2/gal2_length,pgm1=pis_pgm1/pgm1_length,gal1_10_7=pis_gal1_10_7/7296)

p1 = pis_df %>% ggplot(aes(x=gal1_10_7)) + geom_histogram() + xlim(c(0,0.03)) + geom_vline(xintercept = pi_gal/7296, col="red") + theme_bw() + xlab("Nucleotide diversity") +ylab("Count") +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))

p2 = pis_df %>% ggplot(aes(x=pgm1)) + geom_histogram() + xlim(c(0,0.03)) + geom_vline(xintercept = pi_pgm1/pgm1_length, col="red") + theme_bw() + xlab("Nucleotide diversity") +ylab("Count") +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))

p3 = pis_df %>% ggplot(aes(x=gal2)) + geom_histogram() + xlim(c(0,0.03)) + geom_vline(xintercept = pi_gal2/gal2_length, col="red") + theme_bw() + xlab("Nucleotide diversity") +ylab("Count") +
  theme(axis.text = element_text(size=15)) + theme(axis.title= element_text(size=15))

svg("figures/S12.svg",width=12, height=6)
gridExtra::grid.arrange(p1,p2,p3,ncol=3)
dev.off()
