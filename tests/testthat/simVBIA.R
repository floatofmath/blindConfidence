## Test simVBIA function - test more and harder!
context('simVBIA')

test_that("Simualtion of various biases",{
    skip_on_cran()
    ## mean bias under the null is null
    expect_equivalent(round(simVBIA(0,1,1,10^6,runs=10^5)['mean.bias'],4),0)
    ## variance bias for very large effects is null
    expect_equivalent(round(simVBIA(100000,1,1,10,runs=10^5)['variance.bias'],4),0)
    ## variance bias negativen
    expect_less_than(simVBIA(0,1,1,10,runs=10^5)['variance.bias'],0)
    ## negative mean bias positive
    expect_more_than(simVBIA(-1,1,1,10,runs=10^5)['mean.bias'],0)
    expect_less_than(simVBIA(1,1,1,10,runs=10^5)['mean.bias'],0)
    ## absolute unadjusted bias smaller than adjusted
    run <- simVBIA(.5,1,1,8,runs=10^5)
    expect_less_than(run['mean.bias'],run['uc.mean.bias'])
    expect_less_than(run['variance.bias'],run['uc.variance.bias'])
})
