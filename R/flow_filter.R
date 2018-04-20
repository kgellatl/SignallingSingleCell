#' Construct Expression Set Class
#'
#' This function will take an input expression matrix and make an Expression Set Class
#'
#' @param input the input data ex_sc
#' @param panels the list of panels to be gated on
#' @param title the title of the plot
#' @export
#' @details
#' This will take an expression set class and a list of flow like panels and assign cells to groups based on the expression of these genes
#' @examples
#' flow_filter(input = sc_dat, panels = panel_list)

flow_filter <- function(input, panels, title){
  cells <- list()
  for(i in 1:length(panels)){
    panel_n <- panels[[i]]
    tmpdat <- exprs(input)[panel_n,]
    tmpdat[tmpdat > 0] <- 1
    sums <- apply(tmpdat, 2, sum)
    pass_gate <- which(sums == length(panel_n))
    cells_passgate_n <- colnames(tmpdat)[pass_gate]
    cells[[i]] <- cells_passgate_n
  }
  cells_tab <- table(unlist(cells))
  multi_pos <- names(which(cells_tab > 1))
  for(i in 1:length(cells)){
    celltmp <- cells[[i]]
    remove <- match(multi_pos, celltmp)
    remove <- remove[!is.na(remove)]
    if(length(remove > 0)){
      celltmp2 <- celltmp[-remove]
      cells[[i]] <- celltmp2
    }
  }
  pData(input)$Pass_Gate <- NA
  for(i in 1:length(cells)){
    pass_gate_n <- cells[[i]]
    index <- match(pass_gate_n, rownames(pData(input)))
    pData(input)$Pass_Gate[index] <- names(panels)[i]
  }
  pData(input)$Flow_Panel <- NA
  for(i in 1:length(panels)){
    panel <- panels[i]
    ind <- which(pData(input)$Pass_Gate == names(panel))
    panel <- panels[[i]]
    pData(input)$Flow_Panel[ind] <- paste0("C", i, "-", panel, collapse = "\n")
  }
  dat <- pData(input)$Pass_Gate
  gray_points <- which(is.na(dat) == TRUE)
  g <- ggplot()
  g <- g +  theme_classic()
  g <- g +  theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g +  ggtitle(title)
  g <- g +  labs(col= "Panel", x = "tSNE[1]", y = "tSNE[2]")
  g <- g +  geom_point(data = pData(input)[gray_points,c("x", "y")], aes( x=x, y=y), col = "gray", alpha = 0.25)
  g <- g +  geom_point(data = pData(input)[-gray_points,c("x", "y", "Flow_Panel")], aes(x=x, y=y, col=Flow_Panel), alpha = 0.85)
  plot(g)
  return(input)
}

