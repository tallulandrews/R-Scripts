
#library("biomaRt")
#hensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#mensembl = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")


#attr = listAttributes(hensembl)
#hgene_ensemblID = 1
#mouse_orthoID = 848
#hgene_symbol = 63
#entrezID = 57

#filt = listFilters(hensembl)
#biotype = 188

#ensg_name_map <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters = "biotype", values="protein_coding", mart=hensembl)
#ensg2musg <- getBM(attributes=c("ensembl_gene_id","mmusculus_homolog_ensembl_gene"),filters = "biotype", values="protein_coding", mart=hensembl)
#musg_name_map <- getBM(attributes=c("ensembl_gene_id","mgi_symbol"), filters = "biotype", values="protein_coding", mart=mensembl)

#saveRDS(list(ensg_name_map=ensg_name_map, ensg2musg=ensg2musg, musg_name_map=musg_name_map), file="/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/ID_Maps.rds")

obj <- readRDS("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/ID_Maps.rds")
ensg_name_map <- obj[["ensg_name_map"]]
ensg2musg <- obj[["ensg2musg"]]
musg_name_map <- obj[["musg_name_map"]]

h_dups = unique(ensg2musg[(duplicated(ensg2musg[,1])),1])
m_dups = unique(ensg2musg[(duplicated(ensg2musg[,2])),2])
ensg2musg_one2one <- ensg2musg[!(ensg2musg[,1] %in% h_dups) & !(ensg2musg[,2] %in% m_dups),]

map_symbol_ensg <- function(genes, is.org=c("Hsap","Mmus"), is.name=c("symbol","ensg")) {
	if (is.org[1] == "Hsap") {
		if(is.name[1]=="symbol") {
			new = as.character(ensg_name_map[match(genes, ensg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(ensg_name_map[match(genes, ensg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else if (is.org[1] == "Mmus") {
		if(is.name[1]=="symbol") {
			new = as.character(musg_name_map[match(genes, musg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(musg_name_map[match(genes, musg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else {
		stop("Unrecognized organism");
	}
#	new[is.na(new)] = as.character(genes[is.na(new)])
	new[is.na(new)] = ""
	return(new);
}

map_Hsap_Mmus <- function(genes, is.org=c("Hsap","Mmus"), one2one=FALSE) {
	if (one2one) {
		local_map <- ensg2musg_one2one
	} else {
		local_map <- ensg2musg
	}
	if (is.org[1] == "Hsap") {
		new = as.character(local_map[match(genes, local_map[,1]),2])
	} else if (is.org[1] == "Mmus") {
		new = as.character(local_map[match(genes, local_map[,2]),1])
	} else {
		stop("Unrecognized organism");
	}
#	new[is.na(new)] = as.character(genes[is.na(new)])
	new[is.na(new)] = ""
	return(new);
}

General_Map <- function(genes, in.org=c("Hsap","Mmus"), in.name=c("symbol","ensg"), out.org=c("Hsap","Mmus"), out.name=c("symbol","ensg")) {
	if (in.org == out.org & in.name == out.name) {
		# No change
		return(genes)
	}
	if (in.org == out.org) {
		# change names not spp
		return(map_symbol_ensg(genes, is.org=in.org, is.name=in.name))

	} else {
		if (in.name == "symbol") {
			tmp <- map_symbol_ensg(genes, is.org=in.org, is.name=in.name)
		} else {
			tmp <- genes
		}
		tmp <- map_Hsap_Mmus(tmp, is.org=in.org)
		if (out.name =="symbol") {
			out <- map_symbol_ensg(tmp, is.org=out.org, is.name="ensg")
		} else {
			out <- tmp
		}
		return(out)	
	}
}

M_symbol_2_H_symbol <- function(genes) {
        tmp <- map_symbol_ensg(genes, is.org="Mmus", is.name="symbol")
        tmp <- map_Hsap_Mmus(tmp, is.org="Mmus")
        tmp <- map_symbol_ensg(tmp, is.org="Hsap", is.name="ensg")
        return(tmp)
}

H_symbol_2_M_symbol <- function(genes) {
        tmp <- map_symbol_ensg(genes, is.org="Hsap", is.name="symbol")
        tmp <- map_Hsap_Mmus(tmp, is.org="Hsap")
        tmp <- map_symbol_ensg(tmp, is.org="Mmus", is.name="ensg")
        return(tmp)
}

