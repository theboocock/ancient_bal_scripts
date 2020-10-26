get_gff_for_plotting_dn_ds = function(dn_ds,gff)
{
  gff_all = merge(gal1_10_7_dn,gff,by.x="V1",by.y="ID")
  gff_all_sort = gff_all[order(gff_all$order),]
  gff_all_sort$order = 1:nrow(gff_all_sort)
  gff_all_sort$V1 = factor(gff_all_sort$V1)
  
  create_indexes = lapply(split(gff_all_sort, f=gff_all_sort$V1),function(x){x$pos = x$start + seq(from=0,to=length(x$start)*30-1,by=30);return(x)})
  gff_all_sorted_idx = do.call("rbind",create_indexes)
}

get_lod_from_r = function(r,n=365){
  return(-n * (log(1- r^2) / (2 *log(10))))
}

get_gff_for_plotting_dn_ds.2 = function(dn_ds,gff)
{
  gff_all = merge(gal1_10_7_dn,gff,by.x="V1",by.y="ID")
  gff_all_sort = gff_all[order(gff_all$order),]
  gff_all_sort$order = 1:nrow(gff_all_sort)
  gff_all_sort$V1 = factor(gff_all_sort$V1)
  
  create_indexes = lapply(split(gff_all_sort, f=gff_all_sort$V1),function(x){x$pos = x$start + seq(from=0,to=length(x$start)*30-1,by=30);return(x)})
  gff_all_sorted_idx = do.call("rbind",create_indexes)
}