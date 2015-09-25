library(testthat)
library(blindConfidence)

expect_rmonotone <- function(x,y,data,decreasing,digits=4){
    expect_equivalent(round(cor(data[[x]][order(data[[y]])],1:length(data[[x]]),method='spearman'),digits),ifelse(decreasing,-1,1))
}

expect_rincreasing <- function(x,y,data,digits=4){
    expect_rmonotone(x,y,data,decreasing=FALSE,digits=digits)
}
expect_rdecreasing <- function(x,y,data,digits=4){
    expect_rmonotone(x,y,data,decreasing=TRUE,digits=digits)
}


test_check("blindConfidence")
