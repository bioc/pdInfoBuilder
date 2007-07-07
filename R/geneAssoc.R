
fields = c("tx_acc", "role", "distanceToGene", "UnigeneClustID",
   "geneName", "entrezID", "genbankDesc")

clblank = function(x) gsub("^ | $", "", x)

doAssoc = function(x) {
 tmp = strsplit(x, "///")
 lapply(tmp, function(x) { tt = strsplit(x, "//");
         lapply(tt, function(z) {names(z) = fields; gsub("^ | $", "",z)})})
}
