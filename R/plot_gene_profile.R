#' Visualize Profile for a Given Gene
#'
#' @param gene_symbol The symbol of the gene whose profile is to be visualized (character)
#' @param TTCA_result Result of TTCA analysis (data.frame)
#' @param grp1        Data set with longitudinal sampled data (data.frame)
#' @param grp2        Data set with longitudinal sampled data for comparison (data.frame)
#' @param grp1.time   Time points for data set 1 (vector like: c(0,0,0.5,1,2,4,6,8,12,12)
#' @param grp2.time   Time points for data set 2 (vector like: c(0,0,0.5,3,2,4,6,8,12,12,24)
#' @param annot       Annotation for pictures and result (Data.frame with 2 columns with ID and GeneName). (Default: annot=NA)

#' @return displays a plot of the gene expression profile for the gene. Returns the data frame TTCA p-values for the gene (if it exists, if not, returns
#' "Gene not analyzed by TTCA")
#'
#'
#' @importFrom graphics axis
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics title
#'
#'
#' @import MASS
#' @import Matrix
#' @import RISmed
#' @import quantreg
#' @export
visualize_profile <-  function (gene_symbol, TTCA_result, grp1, grp1.time, grp2, grp2.time, annot = NA) 
{
    grp1 <- grp1[, order(grp1.time)]
    grp1.time <- grp1.time[order(grp1.time)]
    grp2 <- grp2[, order(grp2.time)]
    grp2.time <- grp2.time[order(grp2.time)]
    
    
    beginSize = nrow(grp1)
    grp1.grp2 <- merge(grp1, grp2, by = "row.names")
    rownames(grp1.grp2) <- grp1.grp2[, "Row.names"]
    grp1.grp2 <- grp1.grp2[, -1]
    grp1.grp2.time <- c(grp1.time, grp2.time)
    grp <- c(rep(1, length(grp1.time)), rep(2, length(grp2.time)))
    
    rownames(annot) = annot[, "probeset_id"]
    NotAnnot <- rownames((annot[is.na(annot[, "gene_name"]), 
    ]))
    grp1 <- grp1[!(rownames(grp1) %in% NotAnnot), ]
    grp2 <- grp2[!(rownames(grp2) %in% NotAnnot), ]
    grp1.grp2 <- grp1.grp2[!(rownames(grp1.grp2) %in% NotAnnot), 
    ]
    annot <- annot[!(rownames(annot) %in% NotAnnot), ]
    print(paste0(format(Sys.time(), "%H:%M:%S"), ": ", beginSize - 
                     dim(grp1)[1], " not annotated features have been removed (absolute ", 
                 round((1 - dim(grp1)[1]/beginSize) * 100, 1), " % reduction)"))
    
    info <- nrow(grp1)
    Z1 <- apply(grp1.grp2, 1, max)
    if (min(Z1, na.rm = TRUE) < 0) {
        Z2 <- names(Z1[Z1 < 0])
        grp1 <- grp1[!(rownames(grp1) %in% Z2), ]
        grp2 <- grp2[!(rownames(grp2) %in% Z2), ]
        grp1.grp2 <- grp1.grp2[!(rownames(grp1.grp2) %in% Z2), 
        ]
        print(paste0(format(Sys.time(), "%H:%M:%S"), ": ", info - 
                         dim(grp1)[1], " features have been removed because they have no positive intensity values (absolute ", 
                     round((1 - dim(grp1)[1]/beginSize) * 100, 1), " % reduction)"))
        rm(Z2)
    }
    rm(Z1, info)
    
    ## Pot profile
    name <- gene_symbol
    if (nchar(name) > 40) {
        cex1 = 40/nchar(name) * 1.2
    } else {
        cex1 = 1.2
    }
    
    if (!gene_symbol %in% annot$gene_name)
        stop("The gene ", dQuote(gene_symbol), " not in the expression data")
    
    probe_id <- annot$probeset_id[annot$gene_name == gene_symbol][1]
    
    pos <- which(rownames(grp1) == probe_id)
    Poi1 <- as.numeric(grp1[pos, ])
    Poi2 <- as.numeric(grp2[pos, ])
    if (is.na(sum(Poi1, Poi2))) 
        return("at least one of grp1/grp2 contains missing value for the chosen gene")
    
    tmp_df <- data.frame(Poi1 = as.numeric(Poi1),
                         grp1.time = grp1.time)
    p1 <- rqss(Poi1 ~ qss(grp1.time, lambda = 0.6), tau = 0.5, data = tmp_df)
    p1$coef[2:length(p1$coef)] = p1$coef[2:length(p1$coef)] +  p1$coef[1]
    p1 <- p1$coef
    
    tmp_df <- data.frame(Poi2 = as.numeric(Poi2),
                         grp2.time = grp2.time)
    p2 <- rqss(Poi2 ~ qss(grp2.time, lambda = 0.6), tau = 0.5, data = tmp_df)
    p2$coef[2:length(p2$coef)] = p2$coef[2:length(p2$coef)] +  p2$coef[1]
    p2 <- p2$coef
    gl <- c(Poi1, Poi2, p1, p2)
    
    
    gl2 <- max(gl, na.rm = TRUE) - min(gl, na.rm = TRUE)
    
    plot(unique(grp1.time), p1, ylim = c(min(gl) - 
                                             (0.1 * gl2), max(gl, na.rm = TRUE)), type = "l", 
         col = "red", lwd = 2, ylab = "Intensity", xlab = " time [h]", 
         xaxt = "n")
    
    points(grp1.time, Poi1, col = "red")
    lines(unique(grp2.time), p2, col = "blue", lwd = 2)
    points(grp2.time, Poi2, col = "blue", pch = 4)
    axis(side = 1, at = round(unique(grp1.time, grp2.time)))
    title(main = name, cex = cex1)
    legend("bottom", "groups", c("case", 
                                 "control"), pch = c(1, 4), 
           col = c("red", "blue"), ncol = 2, 
           bty = "n")
    
    idx <- which(TTCA_result$gene_name == gene_symbol)
    if (length(idx) == 0)
        return("Gene not analyzed by TTCA")
    
    signif_df <- data.frame(MaxDist = TTCA_result$Pval_MaxDist[idx],
                            Dynamic = TTCA_result$Pval_Dynamic[idx],
                            Early = TTCA_result$Pval_earlyResponse[idx],
                            Middle = TTCA_result$Pval_middleResponse[idx],
                            Late = TTCA_result$Pval_lateResponse[idx],
                            Complete = TTCA_result$Pval_completeResponse[idx],
                            Consensus = TTCA_result$Pval_ConsensusScore[idx])
    signif_df <- as.data.frame(t(signif_df))
    colnames(signif_df) <- gene_symbol
    return(signif_df)
}
    