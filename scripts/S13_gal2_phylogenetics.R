library(ggtree)
library(phangorn)
library(treeio)
library(tidyverse)
sliding_divergence = read.table("data/popgen_in/sliding_divergence_noah.txt", header=F)
#sliding_divergence %>% ggplot(aes(y=1-V3,x=V4 ,color=V1)) + geom_line()
gal2_prot = read.table("data/popgen_in/gal2_protein_structure.txt", header=F,stringsAsFactors = F)
gal2_prot$domain = gal2_prot$V3
gal2_prot$start = gal2_prot$V4 * 3 -2

gal2_prot$end = gal2_prot$V5 * 3 -2
gal2_prot$domain[1] = "N-terminal Cytosolic"
gal2_prot$domain[gal2_prot$domain == "TMhelix"]  = "Transmembrane helix"

gal2_prot$gene = "GAL2"

gal2_prot$domain[gal2_prot$domain == "inside"] = "Cytosolic"
gal2_prot$domain[gal2_prot$domain == "outside"] = "Extracellular"

gal2_prot$domain[nrow(gal2_prot)] = "C-terminal Cytosolic"

for(i in 1:nrow((gal2_prot))){
  if( i==1){
    gal2_prot$end[i] = gal2_prot$end[i] +3
  }else{
    gal2_prot$start[i] = gal2_prot$start[i] + 3
    gal2_prot$end[i] = gal2_prot$end[i]  + 3
  }
}

gal2_prot_all = gal2_prot[1,]

gal2_prot_all$end = tail(gal2_prot$end,n=1)
sliding_divergence$V1 = factor(sliding_divergence$V1,levels=c("gal2a_cbs","gal2a_chin","Suva_10.164","gal2a_baya"))
#sliding_divergence  %>% filter(V2 == "MISSING") %>% ggplot(aes(y=1-V2,x=V3 ,color=V1)) + geom_line() + theme_bw()
#p1 = sliding_divergence  %>% filter(V2 == "MISSING") %>% ggplot(aes(y=abs(V3 * 100),x=(V4 * 10)/3 ,color=V1)) + geom_line() + theme_bw() + ylab("Sequence similarity (%)") + geom_vline(xintercept = 68) +
#  xlab("Base position") + scale_color_manual(values = c("#e41a1c","#377eb8","#984ea3","#ff7f00")) + theme(legend.position = "none") + xlab("Position in gene") + coord_cartesian(xlim=c(0,566.66))
p1 = sliding_divergence  %>% filter(V2 == "MISSING") %>% ggplot(aes(y=abs(V3 * 100),x=(V4 * 10)/3 ,color=V1)) + geom_line() + theme_bw() + ylab("Sequence similarity (%)") + geom_vline(xintercept = 68) +
  xlab("Base position") + scale_color_manual(values = c("#e41a1c","#377eb8","#984ea3","#ff7f00"))  + xlab("Position in gene") + coord_cartesian(xlim=c(0,566.66))

#ggplot(data=gal2_prot,aes(y=))

#ggplot(aes(y=))

gal2_prot$start_gene = gal2_prot_all$start
gal2_prot$end_gene = gal2_prot_all$end
gal2_prot$region = gal2_prot$domain
gal2_prot$region[1] = "N-terminal Cytosolic"
gal2_prot$region[2] = "Rest of the protein"
gal2_prot$end = tail(gal2_prot$end,n=1)
gal2_prot = head(gal2_prot,n=2)
p3 = ggplot() + geom_subgene_arrow(arrowhead_height = unit(12,"mm"),arrowhead_width = unit(12,"mm"),arrow_body_height = unit(8,"mm"),data=gal2_prot,aes(xmin=start_gene/3,xmax=end_gene/3,xsubmin=start/3,xsubmax=end/3,y=gene,fill=domain)) +
  theme_bw() + theme(legend.position = "bottom") + xlab("Position in gene") + coord_cartesian(xlim=c(0,566.66)) + scale_fill_manual(values=c("#4daf4a","grey50"))

svg("figures/S13a.svg",width=5,height=3)
cowplot::plot_grid(p1,p3,nrow=2,rel_heights = c(1,1),axis="lr",align="v")
dev.off()

#strings = readAAStringSet("popgen_phylo_notebooks/data/ygob_w_cbs/R_codon_alignment.fasta")



#ggplot() + geom_subgene_arrow(data=gal2_prot,aes(xmin=start_gene/3,xmax=end_gene/3,xsubmin=start/3,xsubmax=end/3,y=gene,fill=domain)) +
 # theme_bw() + theme(legend.position = "bottom") + xlab("Position in gene") + coord_cartesian(xlim=c(0,566.66))




thatdna_align  = read.raxml("data/trees//RAxML_bipartitionsBranchLabels.R_codon_alignment.fasta.raxml")
p1 = ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("DNA (amino acids 1-67)")
dna_align  = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.R_codon_alignment.aa.raxml")
p2 = ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("Protein (amino acids 1-67)")
dna_align  = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.third.fasta.raxml")
p3=  ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("DNA (amino acids 1-67) 3-fold sites")

dna_align  = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.R_codon_alignment_end.fasta.raxml")
p4 = ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("DNA (amino acids 68-575)")
dna_align  = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.R_codon_alignment_end.aa.raxml")
p5 =ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("Protein (amino acids 68-575)")
dna_align  = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.R_codon_alignment_full.aa.raxml")
p6= ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("DNA full length")
dna_align  = read.raxml("data/trees/RAxML_bipartitionsBranchLabels.R_codon_alignment_full.fasta.raxml")
p7= ggtree(dna_align) +   geom_tiplab()+  geom_nodelab(aes(label=bootstrap,x=branch),vjust=-0.5)  + theme_tree2() + ggtitle("Protein full length")

plots_dna = list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7)
cowplot::save_plot("figures/S13bc.svg",cowplot::plot_grid(p1,p4,nrow=1,rel_heights = c(1,1), align="v"), base_height = 5)

