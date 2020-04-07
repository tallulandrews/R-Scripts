require("biomaRt")
ensembl = useMart("ensembl", mirror='useast')
ensg = useDataset("hsapiens_gene_ensembl", ensembl)
musg = useDataset("mmusculus_gene_ensembl", ensembl)

attributes = listAttributes(musg)


ensg2musg <- getBM(attributes=c('ensembl_gene_id', "mmusculus_homolog_ensembl_gene"), mart = ensg)
ensg2musg <- ensg2musg[ensg2musg[,2] != "",]

ensg2symbol <- getBM(attributes=c('ensembl_gene_id', "hgnc_symbol"), mart = ensg)
ensg2symbol <- ensg2symbol[ensg2symbol[,2] != "",]


musg2symbol <- getBM(attributes=c('ensembl_gene_id', "mgi_symbol"), mart = musg)
musg2symbol <- musg2symbol[musg2symbol[,2] != "",]

saveRDS(list(ensg2musg=ensg2musg, ensg2symbol=ensg2symbol, musg2symbol=musg2symbol), "my_ensembl_name_maps.rds")