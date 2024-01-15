
setwd("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data")

f<-function(count){
  paste0(paste0(which(count!=0),"_",count[which(count!=0)]),collapse="#")
}

list.files("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data",pattern = 'rData') -> file_name

for (i in file_name) {
  load(paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/",i,collapse = ""))
  strsplit(i,split = ".",fixed = T)[[1]][1] -> file_label
  paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/small_size/",file_label,".txt") -> outfile_name
  apply(count,1,f) -> index_data
  rbind(index_data,file_label) -> index_data
  write.table(t(index_data),outfile_name,quote = F)
  print(paste0(file_label," has Finish"))
  gc()
}


for (i in file_name) {
  load(paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/",i,collapse = ""))
  strsplit(i,split = ".",fixed = T)[[1]][1] -> file_label
  print(paste0(dim(count)[2]," ",file_label))
  rm(count)
  gc()
}


list.files("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data",pattern = 'txt') -> file_name

index_count = 0
for (i in file_name) {
  index_count = index_count + 1
  index_data_coord <- read.delim(paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/",i,collapse = ""))
  
  strsplit(file_name[index_count] ,split = ".",fixed = T)[[1]][1] -> file_label
  paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/sup_data/",file_label,"_coord.txt") -> outfile_name
  
  index_data_coord$data_name = file_label
  write.table(index_data_coord,outfile_name,quote = F,row.names = F,sep = "\t")
  print(paste0(file_label," has Finish"))
  gc()
}


index_count = 0
for (i in file_name) {
  index_count = index_count + 1
  index_data_coord <- read.delim(paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/",i,collapse = ""))
  
  strsplit(file_name[index_count] ,split = ".",fixed = T)[[1]][1] -> file_label
  paste0("I:/2021_UNFINSH/exmdb2022/GSE_data/scData/AA_exmdb_data/table_data/",file_label,"_table.txt") -> outfile_name
  
  table(index_data_coord[,c("Cell_type","Sample")])-> a
  as.data.frame.array(a)->data_mix
  apply(data_mix, 2, function(x)(x/sum(x))) -> data_mix
  round(data_mix,4)->data_mix
  row.names(data_mix) -> mix
  gsub("+","_",mix,fixed = T) -> mix
  cbind(mix,data_mix) -> data_mix
  write.table(data_mix,outfile_name,quote = F,row.names = F,sep = "\t")
  print(paste0(file_label," has Finish"))
  gc()
}

data.frame(do.call(rbind,strsplit(rownames(gene_percent_data2),": "))) -> percent_gene_dataset
# split(percent_gene_dataset, list(percent_gene_dataset$X2)) -> split_data_lable
percent_gene_dataset$new_label <- paste0(percent_gene_dataset$X2,"_count")

apply(as.matrix(gene_percent_data2), 1, f) -> scRNA_data_res


rbind(scRNA_data_res,percent_gene_dataset$X1,percent_gene_dataset$new_label) -> cell_exp
t(cell_exp)->cell_exp
as.data.frame(cell_exp) -> cell_exp
rownames(cell_exp)<-NULL
