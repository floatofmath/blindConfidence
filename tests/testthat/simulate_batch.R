context("simulate_batch")
options(mc.cores = min(10,detectCores()-2))
## comparing different applies is meaningless as setting the seed in mclapply is ineffective unless you also change the random number generator
## set.seed(12345)
## mc <- simulate_batch(G,10^6)
## set.seed(12345)
## mc2 <- simulate_batch(G,10^6,TRUE)
## set.seed(12345)
## sc <- simulate_batch(G,10^6,multicore=FALSE)
## mc3 <- simulate_batch(G,10^6)
## set.seed(12345)

test_that("Simulation across a grid of parameters",{
    skip_on_cran()
    G <- expand.grid(delta = c(-1,0,1),sigma = 1,d = 1, n1 = c(2,10,10^6))
    mc <- simulate_batch(G,10^5)
    expect_rincreasing('mean.bias','n1',mc[delta == 1,])
    expect_rdecreasing('mean.bias','n1',mc[delta == -1,])
})





