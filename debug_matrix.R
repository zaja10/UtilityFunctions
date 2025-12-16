library(Matrix)
library(UtilityFunctions)

# Reproduce the matrix
G <- matrix(c(1.00, 0.98, 0.98, 1.00), nrow = 2)
s_obj <- as(Matrix::forceSymmetric(G), "dsCMatrix")

cat("Class of s_obj:\n")
print(class(s_obj))

cat("\nSummary of s_obj:\n")
summ <- summary(s_obj)
print(class(summ))
print(summ)

cat("\nAttempting to access summ$i:\n")
tryCatch(
    {
        print(head(summ$i))
    },
    error = function(e) print(e)
)


cat("\nAttempting conversion to TsparseMatrix:\n")
t_obj <- as(s_obj, "TsparseMatrix")
cat("Class of t_obj:\n")
print(class(t_obj))
cat("Slots:\n")
print(str(t_obj))

cat("\nAccessing slots:\n")
print(t_obj@i + 1)
print(t_obj@j + 1)
print(t_obj@x)
