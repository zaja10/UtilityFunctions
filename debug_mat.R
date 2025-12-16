devtools::load_all()
library(Matrix)

mat <- matrix(c(
    1.0, 0.5, 0.2,
    0.5, 1.0, 0.1,
    0.2, 0.1, 1.0
), nrow = 3, byrow = TRUE)

print("Full Matrix:")
print(mat)
res <- mat2sparse(mat)
print("Sparse Result:")
print(res)

# Debug specific case
mat_sparse <- matrix(c(1, 0, 0, 1), 2, 2)
mat_sparse[upper.tri(mat_sparse)] <- 0
m <- Matrix::Matrix(mat_sparse, sparse = TRUE)
print("Matrix object:")
print(m)
print(class(m))

dsT <- as(as(m, "dgCMatrix"), "TsparseMatrix")
print("dsT object:")
print(dsT)
print(class(dsT))
print("Slots:")
print(dsT@i)
print(dsT@j)
print(dsT@x)

out <- data.frame(
    Row = dsT@i + 1,
    Col = dsT@j + 1,
    Value = dsT@x
)
print("Pre-filter Dataframe:")
print(out)

out <- out[out$Row >= out$Col, ]
print("Post-filter Dataframe:")
print(out)
