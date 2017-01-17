# Tests for MethylMix Package

data(METcancer)
data(METnormal)
data(GEcancer)

results <- MethylMix(METcancer, GEcancer, METnormal)

test_that("Several tests", {
    
    expect_error(MethylMix(METcancer = METcancer, METnormal = METnormal), 
                 'Need to provide GEcancer matrix', fixed = TRUE)
    
    expect_error(MethylMix(GEcancer = GEcancer, METnormal = METnormal), 
                 'Need to provide METcancer matrix', fixed = TRUE)
})
