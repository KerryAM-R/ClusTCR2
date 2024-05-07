
require(ClusTCR2)
example <- read.csv(system.file("extdata", "my_data.csv",package = "ClusTCR2"))

step1 <- ClusTCR(example,allele = F)

step2 <- mcl_cluster(step1)

netplot_ClusTCR2(step2)
