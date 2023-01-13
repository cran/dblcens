UNITEST_d001ch <- function(){
    result = d011ch(z=c(1,2,3,4,5), d=c(1,0,2,2,1), K=3.5, konst=0.6)

    return(sprintf('%0.5f', result[["jump"]][1]));

}

UNITEST_d001 <- function(){
    result = d011(z=c(1,2,3,4,5), d=c(1,0,2,2,1))

    return(sprintf('%0.5f', result[["jump2"]][4]));

}


test_that("dblcens works", {
  expect_equal(UNITEST_d001ch(), "0.49996")
  
  expect_equal(UNITEST_d001(), "0.29999")
})
