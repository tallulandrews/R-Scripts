Complex_Markers <- function (expr_mat, labels, n_max=1, strict_only = FALSE) {
        # n_max = 1 should give the same as original just using different package.
        require("pROC")
        if (length(labels) != length(expr_mat[1, ])) {
                stop("Length of labels does not match number of cells.");
        }
        if (n_max > (length(unique(labels))-1)) {
                stop("n_max must be less than the number of labels");
        }
	if (min(factor_counts(labels)) < 10) {
		print("Warning: Small groups (n < 10) may bias marker gene results.")
	}
	# Remove groups with a single cell
	label_counts = summary(factor(labels))
	exclude = names(label_counts)[label_counts<2]
	if (length(exclude) > 0) {
		warning(paste("Warning: Excluding",length(exclude),"groups with less than 2 samples."))
#		print(paste("Warning: Excluding",length(exclude),"groups with less than 2 samples."))
		keep = !(labels %in% exclude)
		expr_mat <- expr_mat[,keep]
		labels <- labels[keep]
		labels <- factor(labels);
	}

        # Mean ranked expression all genes, each cluster (efficient)
        gene_cluster_means <- function(mat, groups) {
                MAT <- as.matrix(mat)
                x <- split(seq(ncol(MAT)), groups)
                result <- sapply(x, function(a) rowMeans(MAT[,a]))
                return(result);
        }
        ranked_matrix <- t(apply(expr_mat, 1, rank))
        gene_cluster_ranks <- gene_cluster_means(ranked_matrix, labels)
        cluster_priority <- t(apply(-1*gene_cluster_ranks, 1, rank))
	label_columns <- sort(unique(labels));
        # Consider each gene a possible marker for top 1->n_max groups
        gene_auc <- function(g) {
		# Get AUC + 95% CI for this gene in top n groups
                get_auc_ci <- function(n) {
                        g_groups <- colnames(cluster_priority)[which(cluster_priority[g,] <= n)]
                        if (length(g_groups) < length(unique(labels)) & length(g_groups) > 0) {
                                g_roc <- roc(controls=expr_mat[g,!(labels %in% g_groups)], cases=expr_mat[g,(labels %in% g_groups)],
                                                auc=TRUE, ci=TRUE, direction="<");
                                return(as.vector(g_roc$ci))
                        } else {
                                return(c(0,0.5,1))
                        }
                }
		# AUCs for this gene at all ns
                auc_tab <- sapply(1:n_max, get_auc_ci)

		# Identify top set of groups this gene is a marker for & determine if good enough to return.
                top = which(auc_tab[2,] == max(auc_tab[2,]))
                sec = which(auc_tab[2,] == max(auc_tab[2,-top]))
                if (n_max > 1 & min(auc_tab[1,top]) <= max(auc_tab[2,sec]) & strict_only) {
                        # Not a marker : ci of top set of groups contains auc of second best set of groups
                        group=rep(0, times=length(label_columns));
                        pval=-1
                        auc=-1
                } else {
                        # Return marker info
                        n = max(top)
                        g_groups = colnames(cluster_priority)[which(cluster_priority[g,] <= n)];
			if (n_max == 1) {
	                        group = paste(g_groups, collapse="+")
			} else {
				group = as.numeric(label_columns %in% g_groups);
			}
                        auc = auc_tab[2,n]
			# p.value from wilcox test
                        pval = wilcox.test(expr_mat[g,!(labels %in% g_groups)],expr_mat[g,(labels %in% g_groups)])$p.value
                }
                return(c(auc, group, pval))
        }
	# Get best AUC across sets of groups for each gene
        out <- sapply(1:length(expr_mat[,1]),gene_auc)

	# Format output nicely
        out_matrix = as.data.frame(t(out))
        rownames(out_matrix) <- rownames(expr_mat)
	if (n_max == 1) {
	        colnames(out_matrix) <- c("AUC", "Group",  "p.value")
	} else {
		colnames(out_matrix) <- c("AUC", as.character(label_columns), "p.values")
	}
        out_matrix[,1] = as.numeric(as.character(out_matrix[,1]))
        out_matrix[,3] = as.numeric(as.character(out_matrix[,3]))
	# Apply bonferroni (-ish) correction
	# considers number of groups & number of genes, but not all possible combinations of groups
        out_matrix$q.value = out_matrix$p.value*length(unique(labels))*length(expr_mat[,1]);
        out_matrix$q.value[out_matrix$q.value < 0] = -1;
        out_matrix$q.value[out_matrix$q.value > 1] = 1;
        return(out_matrix);
}

# Turn matrix of on/off into combined names of positive or negative groups
get_combo_names <- function(marker_matrix) {
	tmp <- marker_matrix[,2:(length(marker_matrix[1,])-2)]
	out <- apply(tmp, 1, function(x){paste(colnames(tmp)[x==1], collapse="+")})
	# If gene is "off" in fewer groups than it is "on" it is a negative marker
	anti = rowSums(tmp) > 0.5*length(tmp[1,])
	out_anti <- apply(tmp, 1, function(x){paste(colnames(tmp)[x==0], collapse="+")})
	out[anti] = paste("NOT", out_anti[anti])
	return(out);
}

factor_counts <- function(labels) {
	labels <- factor(labels)
        x <- split(seq(length(labels)), labels)
        result <- sapply(x, function(a) length(a))
        return(result);
}
marker_heatmap.3 <- function(marker_matrix, expr_mat, cell_colour_bars=NULL, min_AUC=0.7, max_q_value=1, min_cat_size=10) {
	#Get marker categories
	Gene_assign <- get_combo_names(marker_matrix)



	other_cats <- factor_counts(Gene_assign) < min_cat_size
	Gene_assign[Gene_assign %in% names(other_cats[other_cats])] <- "Other"

	# Apply filters
	keep <- rowSums(marker_matrix[,2:(length(marker_matrix[1,])-2)]) > 0 & 
			marker_matrix[,1] > min_AUC & Gene_assign != "Other" &
			marker_matrix[,length(marker_matrix[1,])] <= max_q_value
	toplot <- expr_mat[keep,]
	gene_colour_bar = factor(Gene_assign[keep])
	gene_colour_palette = rainbow(length(levels(gene_colour_bar)));

	# Order genes appropriately
	sort_order <- order(Gene_assign[keep], rowMeans(toplot))
	toplot<- toplot[sort_order,]
	gene_colour_bar <- gene_colour_bar[sort_order]

	# Heatmap Setup
	require("RColorBrewer")
	source("~/R-Scripts/heatmap.3.R")
	heatcolours <- rev(brewer.pal(11, "RdBu"))
	col_breaks <- c(-100, seq(-2, 2, length = 10), 100)
	
	# Heatmap
	heatout <- heatmap.3(toplot, scale="row", breaks=col_breaks, col=heatcolours, Rowv=FALSE,
                        ColSideColors=cell_colour_bars, ColSideColorsSize=length(cell_colour_bars[1,]),
                        RowSideColors=matrix(gene_colour_palette[gene_colour_bar], nrow=1), RowSideColorsSize=1,
                        hclustfun=function(x){hclust(x, method="ward.D2")})

	names(gene_colour_palette) = levels(gene_colour_bar)
	return(list(legend=gene_colour_palette, gene_assignments=Gene_assign));
}
