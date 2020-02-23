load("plots/example1.RData")
methods <- c("Replication ARM", "SOIL ARM", "Replication BIC", "SOIL BIC")
i <- 1
for (method in example1$`Plot  4`) {
  print(methods[i])
  print(which(method > 0.02 ))
  i <- i+1
}

