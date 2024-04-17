#########################################
# wrapper for rule: merge_detected_species
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[[1]]
sink(logfile, append = T, type = "output")
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: merge_detected_species \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(data.table)

intab = snakemake@input[["table"]]
outtab = snakemake@output[['table']]
nreads = as.numeric(snakemake@params[['nreads']])
evalue = as.numeric(snakemake@params[['evalue']])
tmpd = snakemake@params[['tmpd']]
top = snakemake@params[['top_sp']]

tab = data.table()
for(inp in intab) {
  name = gsub('(.+).species_stats.tsv','\\1',basename(inp))
  cat("## Reading detected_species table for sample",name,"into R for processing\n")
  tab = rbind(tab, fread(file = inp, sep = '\t', header = T, col.names = c("Species","#Reads","%Reads"))[,.(Species,`#Reads`,`%Reads`,Sample=name)])
}
cat("## Processing merged table\n")
tab = dcast(tab[,.SD,.SDcols=!c('%Reads')], Species ~ Sample,  fill = 0, value.var = '#Reads')
tab[, max:=max(.SD), by=seq_along(Species)]
tab[, sum:=sum(.SD), by=seq_along(Species)]
tab = tab[order(-sum)][,.SD,.SDcols=!c('max','sum')]

unmapped = c("Species"="Insignificant hits", tab[, nreads-colSums(.SD), .SDcols=!'Species'])
if(tab[,.N] > top) {
  other = c("Species"="Other species", tab[(top+1):.N, colSums(.SD), .SDcols=!'Species'])
  tab = rbind(tab[1:top], as.list(other))
}
tab = rbind(tab, as.list(unmapped))

cat("## Writing summary table down\n")
cat('# id: "Species detection summary table"\n', file = outtab)
cat('# section_name: "Species detection"\n', file = outtab, append = T)
cat('# description: "Species detection against <a href=https://www.ncbi.nlm.nih.gov/nucleotide/>Nucleotide database</a> for randomly subsampled',nreads,'reads with valid hits having E-value <',evalue,'"\n', file = outtab, append = T)
cat('# format: "tsv"\n', file = outtab, append = T)
cat('# plot_type: "bargraph"\n', file = outtab, append = T)
fwrite(as.data.table(t(tab), keep.rownames = T),
       outtab,
       sep = '\t',
       col.names = F,
       row.names = F,
       append = T,
       quote = F)