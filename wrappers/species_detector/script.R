#########################################
# wrapper for rule: species_detector
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[[1]]
sink(logfile, append = T, type = "output")
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: species_detector \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(data.table)

fastq_file = snakemake@input[["fastq"]]
taxdbd = snakemake@input[['taxdbd']]
taxdbi = snakemake@input[['taxdbi']]
outtab = snakemake@output[['table']]
threads= as.numeric(snakemake@threads)
ntdb   = snakemake@params[['ntdb']]
fasta  = snakemake@params[["fasta"]]
blast  = snakemake@params[["blast"]]
nreads = as.numeric(snakemake@params[['nreads']])
evalue = as.numeric(snakemake@params[['evalue']])
tmpd   = snakemake@params[['tmpd']]

cmd = paste0("( time seqtk seq -a ",fastq_file,"| seqtk sample -s 123 - ",nreads," ) > ",fasta," 2>> ",logfile)
cat("## COMMAND: ",cmd,"\n")
shell(cmd)

cmd = paste0("cp ",taxdbd," ./ >> ",logfile," 2>&1")
cat("## COMMAND: ",cmd,"\n")
shell(cmd)

cmd = paste0("cp ",taxdbi," ./ >> ",logfile," 2>&1")
cat("## COMMAND: ",cmd,"\n")
shell(cmd)

cmd = paste0('( time blastn -max_hsps 1 -num_threads ',threads,' -max_target_seqs 1 -db ',ntdb,' -query ',fasta,
             ' -out ',blast,' -outfmt "6 std staxids sscinames scomnames sblastnames" ) >> ',logfile,' 2>&1')
cat("## COMMAND: ",cmd,"\n")
shell(cmd)

# cmd = paste0('( sort -k 1,1 -k 12,12nr $FASTQ_BASE.blast > $FASTQ_BASE.blast.sorted && mv $FASTQ_BASE.blast.sorted $FASTQ_BASE.blast ) 2>> ')
colnames = c("query","target","pident","length","mismatch","gapopen","qstart","qend",
             "tstart","tend","eval","bitscore","ttaxid","Species","tcomname","tblastname")
cat("## Reading BLAST output into R for post-processing\n")
tab = fread(file = blast, sep = '\t', header = F, col.names = colnames)

if(tab[,.N] == 0 || tab[eval <= evalue,.N] == 0) {
  cat("## Writing empty summary table: ",outtab,"\n")
  cat(paste0("# Input file: ",fastq_file,"\n"), file = outtab)
  cat(paste0("# Number of reads with BLAST hits to nt DB (e-value < ",evalue,"): ",0,"/",nreads,"\n"), 
      file = outtab, append = T)
  fwrite(list("Species","#Reads","%Reads"),
        outtab,
        sep = '\t',
        col.names = T,
        row.names = F,
        append = T,
        quote = F)
} else {
  cat("## Writing summary table: ",outtab,"\n")
  cat(paste0("# Input file: ",fastq_file,"\n"), file = outtab)
  cat(paste0("# Number of reads with BLAST hits to nt DB (e-value < ",evalue,"): ",tab[eval <= evalue, .N],"/",nreads,"\n"), 
      file = outtab, append = T)
  fwrite(tab[eval < evalue, .(`#Reads`=.N,`%Reads`=round((.N/nreads)*100, 2)), by=.(Species)][order(-`#Reads`)],
        outtab,
        sep = '\t',
        col.names = T,
        row.names = F,
        append = T,
        quote = F)
}
