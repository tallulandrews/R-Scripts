map = read.table("/lustre/scratch108/compgen/team218/TA/Mirrored_Annotations/Mmus_Gene_Name_Mapping_Ensembl80.out", header=T)
ensg2symbol <- function(x) {
        new = as.character(map[match(x, map[,1]),2])
        new[is.na(new)] = as.character(x[is.na(new)])
        new[duplicated(new)] = x[duplicated(new)]
        return(new)
}
