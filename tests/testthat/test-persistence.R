mat <- matrix(c(0, 1, 2, 1, 0.5, 2, 2, 1, 3), ncol = 3, byrow = TRUE)
pd <- as_persistence(mat)
df <- data.frame(
  dimension = c(1, 1, 1, 1, 1),
  x = c(1, .6, 1.3, 2.1, 1.99),
  y = c(2, 1.2, .7, 4.2, 1.99)
)


test_that("as_persistence.default correctly processes a basic matrix", {
  
  expect_s3_class(pd, "persistence")
  expect_equal(length(pd$pairs), 3)  # 3 unique degrees (0,1,2)
  expect_equal(pd$pairs[[1]], matrix(c(1, 2), ncol = 2))  # Degree 0
  expect_equal(pd$pairs[[2]], matrix(c(0.5, 2), ncol = 2)) # Degree 1
  expect_equal(pd$pairs[[3]], matrix(c(1, 3), ncol = 2))  # Degree 2
})

test_that("as_persistence.default correctly processes output from `____$diag`", {
  coords <- df[, c("x", "y")]
  
  pd <- alphaComplexDiag(df)$diagram
  pd[, c(2, 3)] <- sqrt(pd[, c(2, 3)])
  pd_p <- as_persistence(pd)
  
  df3d <- data.frame(
    x = c(1, 0.6, 1.3, 2.1, 1.99),
    y = c(2, 1.2, 0.7, 4.2, 1.99),
    z = c(0.1, 0.3, 0.2, 0.4, 0.35)
  )
  
  pd2 <- alphaShapeDiag(df3d)$diagram
  pd2[, c(2, 3)] <- sqrt(pd2[, c(2, 3)])
  pd_p2 <- as_persistence(pd2)
  
  FltRips <- ripsFiltration(X = df, maxdimension = 1,
                            maxscale = 1.5, dist = "euclidean", library = "Dionysus",
                            printProgress = FALSE)
  DiagFltRips <- filtrationDiag(filtration = FltRips, maxdimension = 1,
                                library = "Dionysus", location = TRUE, printProgress = TRUE)
  pd3 <- DiagFltRips$diagram
  pd3[, c(2, 3)] <- sqrt(pd3[, c(2, 3)])
  pd_p3 <- as_persistence(pd3)
  
  Diag1 <- gridDiag(coords, distFct, lim = cbind(c(-1, 1), c(-1, 1)), by = 0.05, sublevel = TRUE,
                    printProgress = TRUE) 
  pd4 <- Diag1$diagram 
  pd4[, c(2, 3)] <- sqrt(pd4[, c(2, 3)])
  pd_p4 <- as_persistence(pd4)
  
  pd5 <- ripsDiag(df, maxdimension = 1, maxscale = 10)$diagram
  pd5[, c(2, 3)] <- sqrt(pd5[, c(2, 3)])
  pd_p5 <- as_persistence(pd5)
  
  expect_s3_class(pd_p5, "persistence")
  expect_equal(pd_p5$max_dim, 1)
  expect_equal(length(pd_p5$pairs[[1]][1, ]), 2) # Correct Format
  
  expect_s3_class(pd_p4, "persistence")
  expect_equal(pd_p4$max_dim, 1)
  expect_equal(length(pd_p4$pairs[[1]][1, ]), 2) # Correct Format
  
  expect_s3_class(pd_p3, "persistence")
  expect_equal(pd_p3$max_dim, 1)
  expect_equal(length(pd_p3$pairs[[1]][1, ]), 2) # Correct Format
  
  expect_s3_class(pd_p2, "persistence")
  expect_equal(pd_p2$max_dim, 1)
  expect_equal(length(pd_p2$pairs[[1]][1, ]), 2) # Correct Format
  
  expect_s3_class(pd_p, "persistence")
  expect_equal(pd_p$max_dim, 1)
  expect_equal(length(pd_p$pairs[[1]][1, ]), 2) # Correct Format
  
  
})

test_that("as_persistence.PHom correctly processes output from `vietoris_rips` 
          from {ripserr}", {
            pd <- cubical(volcano)
            pd_p <- as_persistence(pd)
            
            expect_s3_class(pd_p, "persistence")
            expect_equal(pd_p$max_dim, 1)
            expect_equal(length(pd_p$pairs[[2]][1, ]), 2)
          })


test_that("as_persistence.list recongnizes a list containing a matrix", {
  mat <- matrix(c(0, 1, 2, 1, 0.5, 2), ncol = 3, byrow = TRUE)
  pd <- as_persistence(list(mat))
  
  expect_s3_class(pd, "persistence")
  expect_equal(length(pd$pairs), 2)  # Degrees 0 and 1
})


test_that("as_persistence.list handles a list of matrices correctly", {
  pd <- as_persistence(list(matrix(c(1, 2), ncol = 2), matrix(c(0.5, 2), ncol = 2)))
  
  expect_s3_class(pd, "persistence")
  expect_equal(length(pd$pairs), 2)  # Two degrees (0 and 1)
  expect_equal(pd$pairs[[1]], matrix(c(1, 2), ncol = 2))
  expect_equal(pd$pairs[[2]], matrix(c(0.5, 2), ncol = 2))
})


test_that("as_persistence.persistence returns the object unchanged", {
  mat <- matrix(c(0, 1, 2, 1, 0.5, 2), ncol = 3, byrow = TRUE)
  pd <- as_persistence(mat)
  pd2 <- as_persistence(pd)
  
  expect_identical(pd, pd2)
})


test_that("as.data.frame.persistence creates correct format", {
  mat <- matrix(c(0, 1, 2, 1, 0.5, 2), ncol = 3, byrow = TRUE)
  pd <- as_persistence(mat)
  df <- as.data.frame(pd)
  
  expect_s3_class(df, "data.frame")
  expect_equal(colnames(df), c("degree", "birth", "death"))
  expect_equal(as.numeric(df[1, ]), c(0, 1, 2))
})

test_that("`get_pairs` correctly grabs pairs from persistence object", {
  pd_p <- as_persistence(df)
  expect_equal(get_pairs(pd_p, 1), pd_p$pairs[[2]])
  expect_true(all(is.na(get_pairs(pd_p, 3)))) #nonexistent dimension returns empty matrix
  expect_error(get_pairs(pd)) #verifies failure given incorrect class
})

