#'Plotting Significant features data
#'
#'This function takes data generated from differential expression analysis functions (see \link[microbialSeq]{kruskal_expression},
#'\link[microbialSeq]{differential_expression}) and returns ggplot objects showing different visualisation given the available parameters.
#'
#' @param df (optional). A \code{data.frame} returned by \link[microbialSeq]{kruskal_expression},
#'\link[microbialSeq]{differential_expression}
#' @param df_accuracy(optional). A \code{data.frame} of  feature importance information.
#' @param  res_tax (optional). A \code{data.frame} bearing results for plotting MA plot.
#' @param corrections_table (optional). A \code{data.frame} containing multiple testing corrections information.
#' @return Returns ggplot object(s). This can further be manipulated as preferred by user.
#' @examples 
#' 
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' 
#' @author Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#' 
#' @import ggplot2
#' @import grid
#'
#' @export plot_signif
#' 
plot_signif <- function(df=NULL, df_accuracy=NULL, res_tax=NULL,pvalue.Cutoff=0.05, corrections_table=NULL, ...){
  
  #==plot the significant features information and random classifier results 
  p<-NULL
  if(!is.null(df)){
    p<-ggplot(df,aes(Groups,Value,colour=Groups))
    p<-p+geom_boxplot(outlier.size = NA)+geom_jitter(position = position_jitter(height = 0, width=0))+theme_bw()
    p<-p+ facet_grid( ~Rank+Taxa, scales="free_x")
    p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text = element_text(size = 10, colour = "black", angle = 90,vjust=0)) #vjust=0 aligns the labels to the bottom in this case
    p<-p+theme(strip.background =element_rect(fill="white"))+theme(plot.margin = unit(c(1, 1, 0, 1), "lines"))#set background colour for facet lbels
  }
  
  #==generate stand alone plot for random forest results =============#
  prf<-NULL
  if(!is.null(df_accuracy)){
    prf <- ggplot(data = df_accuracy,aes(x=Sample,y=Value)) + theme_bw()
    prf <- prf+geom_bar(stat = "identity",fill="darkblue",width = 0.5)
    prf <- prf + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    prf <- prf + xlab("Taxa description") + ylab("Mean Decrease in Accuracy")
  }
    
  #===produce  MA plot ===============#
  p1 <- NULL
  p2<-NULL
  if(!is.null(res_tax)){
    label=T
    tax.display = NULL
    tax.aggregate = "OTU"
    p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) + geom_point(size = 2) + scale_x_log10()
    p1 <- p1+scale_color_manual(values=c("black", "red"))+labs(x="Mean abundance",y="Log2 fold change")+theme_bw()
    if(label == T){
      if (!is.null(tax.display)){ 
        rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
      } 
      else {
        rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
      }
      p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
    }
  
    #== produce a bar version of ma plot =========================#
    res_tax$cols<- ifelse(res_tax$log2FoldChange>=0, "Up regulated", "Down regulated")
    p2<-ggplot(data=res_tax,aes(x=rownames(res_tax),y=log2FoldChange,fill=cols))+geom_bar(stat="identity",width=0.5)+theme_bw()
    p2<-p2 + theme(axis.text.x = element_text(angle = 90,hjust = 1))+ xlab("Taxa Description") +ylab("Log2 Fold Change")+
      geom_text(aes(label=round(as.numeric(baseMean),1)), vjust=0)
    p2<-p2+scale_fill_manual(values = c("Up regulated" = "darkblue", "Down regulated" = "red"))+theme(legend.title=element_blank())
  }
  
  #== plot of multiple testing corrections ====
  if(!is.null(corrections_table)){
    plot(corrections_table$p.value, corrections_table$E.value,main='Multitesting corrections',
         xlab='Nominal p-value',ylab='Multitesting-corrected statistics',log='xy',col='blue',panel.first=grid(col='#BBBBBB',lty='solid'))
    lines(corrections_table$p.value,corrections_table$FWER,pch=20,col='darkgreen', type='p')
    lines(corrections_table$p.value,corrections_table$q.value,pch='+',col='darkred', type='p')
    abline(h=pvalue.Cutoff, col='red', lwd=2)
    legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
  }
  out <- list("significantfeaturesplot"=p, "MAplot"=p1, "lfcplot"=p2, "MDAplot"=prf)
}
