##' Plot coverage probabilities of confidence intervals following blinded sample size reassessment
##'
##' @title Coverage plot
##' @param gridsim Simulation results
##' @param adjusted should the results for the adjusted interim variance estimate be plotted
##' @return ggplot object
##' @author float
##' @export
plot_coverage <- function(gridsim,adjusted=F){
    if(adjusted){
            tupd <- ggplot(gridsim,mapping=aes(x=delta,y=100*(.025-uc.upper.prob))) +
                geom_line(y=0,colour='gray') +
                geom_path(lty=2) + 
                geom_path(aes(y=100*(.025-uc.lower.prob)),lty=3) + 
                geom_path(aes(y=100*(.05-uc.total.prob)))
    } else {
        tupd <- ggplot(gridsim,mapping=aes(x=delta,y=100*(.025-upper.prob))) +
            geom_line(y=0,colour='gray') +
            geom_path(lty=2) + 
            geom_path(aes(y=100*(.025-lower.prob)),lty=3) + 
            geom_path(aes(y=100*(.05-total.prob)))
    }
    tupd <- tupd + facet_grid(s~sigma,labeller=.lgrid,scales="free_y") + theme_bw() + #ylim(-.025,.01) +
        xlab(expression(delta)) + ylab('Actual - nominal coverage [%]')
    print(tupd)
    return(tupd)
}

##' Plot the bias of the effect estimate following blinded sample size reassessment
##'
##' @title Mean bias plot
##' @param gridsim Simulation results
##' @param thmax_bias Should the theoretical upper bound be plotted (needs to be in the data)
##' @return ggplot object
##' @author float
##' @export
mean_bias_plot <- function(gridsim,thmax_bias=TRUE){
    tupd.mean.bias <- ggplot(gridsim,mapping=aes(x=delta,y=mean.bias)) +
        geom_path(lty=2) + geom_line(y=0,colour='gray') +
        geom_path(aes(y=uc.mean.bias),lty=1) +
        geom_path(aes(y=brannath),lty=2,col='darkgray',data=subset(gridsim,brannath <= 0.1)) +
        geom_path(aes(y=-brannath),lty=2,col='darkgray',data=subset(gridsim,brannath <= 0.1)) +

        facet_grid(s~sigma,labeller=.lgrid) + theme_bw() +
    xlab(expression(delta)) + ylab('Bias of the mean') #+ ylim(-.4,.4)
    if(thmax_bias) {
        tupd.mean.bias <- tupd.meanbias +
                    geom_path(aes(y=m.bias),lty=3,col='black',data=subset(gridsim,m.bias<=0.2 & n1 == 2)) +
            geom_path(aes(y=m.bias.n),lty=3,col='black',data=subset(gridsim,m.bias.n>=-0.2 & n1 == 2)) +
            geom_path(aes(y=m.bias),lty=3,col='black',data=subset(gridsim,m.bias<=0.1 & n1 == 8)) +
            geom_path(aes(y=m.bias.n),lty=3,col='black',data=subset(gridsim,m.bias.n>=-0.1 & n1 == 8)) +
            geom_path(aes(y=m.bias),lty=3,col='black',data=subset(gridsim,m.bias<=0.1 & n1 == 18)) +
            geom_path(aes(y=m.bias.n),lty=3,col='black',data=subset(gridsim,m.bias.n>=-0.1 & n1 == 18)) +
            geom_path(aes(y=m.bias),lty=3,col='black',data=subset(gridsim,m.bias<=0.1 & n1 == 32)) +
            geom_path(aes(y=m.bias.n),lty=3,col='black',data=subset(gridsim,m.bias.n>=-0.1 & n1 == 32))
    }
    print(tupd.mean.bias)
    tupd.mean.bias
}


##' Plot the bias of the variance of the effect estimate following blinded sample size reassessment 
##'
##' @title Variance bias plot
##' @param gridsim Simulation results
##' @return ggplot object
##' @author float
##' @export
variance_bias_plot <- function(gridsim){
    tupd.var.bias <- ggplot(gridsim,mapping=aes(x=delta,y=variance.bias)) +
#    geom_line(y=0,colour='gray',lty=2) +
    geom_path(aes(y=bound),color="darkgray") +
    geom_path(lty=2) + geom_path(aes(y=uc.variance.bias),lty=1) +
    facet_grid(s~sigma,labeller=.lgrid) + theme_bw() + #ylim(-0.3,0.1) +
    xlab(expression(delta)) + ylab('Bias of the variance')
    print(tupd.var.bias)
    tupd.var.bias
}


##' Plot the variance of the effect estimate (the square of the standard error of the mean) following blinded sample size reassessment  
##'
##' @title SEM plot
##' @param gridsim Simulation results
##' @param adjusted should the results for the adjusted interim variance estimate be plotted
##' @return ggplot object
##' @author float
##' @export
adjusted_sem_plot <- function(gridsim,adjusted=F){
    if(adjusted){
        tupd <- ggplot(gridsim,mapping=aes(x=delta,y=uc.vm)) +
            geom_line(lty=1) + geom_path(aes(y=uc.ev),lty=2) +
            geom_path(aes(y=uc.exv),lty=3) 
    } else {
        tupd <- ggplot(gridsim,mapping=aes(x=delta,y=vm)) +
            geom_line(lty=1) + geom_path(aes(y=ev),lty=2) +
            geom_path(aes(y=exv),lty=3)
    }
    tupd <- tupd + theme_bw() +
        xlab(expression(delta)) +
        geom_blank(aes(x=delta,y=ylow)) +
        ylab('Variance of the mean estimate')+
        facet_grid(s~sigma,labeller=.lgrid) 
    print(tupd)
    tupd
}

maximum_bias_plot <- function(maxsim){
    bias_string <- sub("uc\\.","",maxsim$what.max)[1]
    bias_short <- c('mean.bias' = 'str_mmb','variance.bias' = 'str_mvb')[bias_string]
    
    m_vec <- c("adjusted","unadjusted")
    names(m_vec) <- c(bias_string,paste0("uc.",bias_string))
    max.bias2 <- as.data.frame(maxsim)
    max.bias2$method  <- factor(m_vec[maxsim$what.max])
    colnames(max.bias2) <- sub(paste0("(..\\.)?(",bias_string,")"),"\\1max.\\2",colnames(max.bias2))
    max.bias2 <- max.bias2[,c('method','n1','delta','sigma',paste0('max.',bias_string),paste0('uc.max.',bias_string))]
    df <- melt(max.bias2,id=c('n1','method'))
    df <- subset(df,!(method=='adjusted' & variable == paste0('uc.max.',bias_string)))
    df <- subset(df,!(method=='unadjusted' & variable == paste0('max.',bias_string)))
    df <- within(df,scale <- c("absolute","relative")[grepl("rel",variable)+1])
    df <- within(df,variable <- sub("rel\\.","",variable))
    df <- within(df,variable <- sub("uc\\.","",variable))
    df <- df[df$variable %in% c(paste0('max.',bias_string),'delta','sigma'),]
    df <- within(df,variable <- factor(variable))
    rev_vec <- c(bias_short)
    names(rev_vec) <- paste0('max.',bias_string)
    df <- within(df,variable <- revalue(variable,rev_vec))
    df <- within(df,variable <- relevel(variable,'delta'))
    df <- within(df,variable <- relevel(variable,bias_short))
    df <- within(df,method <- relevel(method,'unadjusted'))

    dfa <- df[df$scale=="absolute",]#"relative",]#
    scalehelper <- data.frame(n1=rep(-10,3),variable=factor(1:3,label=c(bias_short,expression(delta),expression(sigma))),value=rep(0,3),scale=rep("absolute",3),method=rep("adjusted",3))
    dfa <- rbind(dfa,scalehelper)
    
    mplot <- ggplot(dfa)+
        geom_line(aes(n1,-value),data=subset(dfa,variable==bias_short)) +
        geom_smooth(aes(n1,value),data=subset(dfa,variable!=bias_short),col='darkgray',fill='gray') + geom_point(aes(n1,value),data=subset(dfa,variable!=bias_short),size=.2) +
        facet_grid(variable~method,scale='free_y',labeller=lverbose) +xlab(expression(n[1]))+ylab('')+theme_bw()+xlim(c(2,50))
    
    print(mplot)
}



.lmax  <- function(labels,from=c("anc","mmb","mvb"),to=c("AC-NC [%]","max. mean bias","abs. variance bias"),multi_line=TRUE){
    str_labels <- label_value(labels,multi_line)
    par_labels <- label_parsed(labels,multi_line)
    strings <- grep("str_",str_labels[[1]])
    str_labels[[1]] <- sub("str_","",str_labels[[1]])
    if(!is.null(from)){
        str_labels[[1]] <- mapvalues(str_labels[[1]],from,to)
    }
    par_labels[[1]][strings] <- str_labels[[1]][strings]
    return(par_labels)
}


.lgrid  <- function(labels,from=c("anc","mmb","mvb"),to=c("actual - nominal coverage [%]","maximum mean bias","absolute bias"),multi_line=TRUE){
    variable <- names(labels)
    my_labels <- list()
    if(variable == "s"){
        my_labels[[1]] <- sapply(labels[[1]],function(val) substitute(paste(sigma[0]," = ",foo,", ",n[1]," = ",n1,sep=''),
                                              list(foo=val,n1=ceiling(1/2*zss(val,1,.025,.2)))))
    } else {
        my_labels[[1]] <- sapply(labels[[1]],function(val) substitute(paste(sigma," = ",foo,sep=""),list(foo=val)))
    }
    return(my_labels)
}
    
