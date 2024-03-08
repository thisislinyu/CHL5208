args = commandArgs(trailing=TRUE)
# You must include the statement above so the "argument" (CHR) is passed into the R script
# Function of this script: Read .raw file and segment into chunks of 25000 each

# Remember: the path that is used to read the data is set by you in the bash script (e.g., MAINDIR)
# This is the input/read in line
geno = read.table(paste0('fracture_chr',args[1],'.qced_recoded_geno.raw'),header=T,
	colClasses='character')

(numvar = ncol(geno) - 6) # number of variants

breaks = seq(1,numvar,50000)
if (numvar>=max(breaks)) { breaks = c(breaks,(numvar+1)) }
print(breaks)

frontpart = geno[,1:6]
str(frontpart)

check = NULL
for (k in 1:(length(breaks)-1)) {
	g = NULL
	start = breaks[k]+6
	end = breaks[k+1]-1+6
	varnames = colnames(geno)[start:end]
	g = as.data.frame(geno[,c(start:end)],drop=F,stringsAsFactors=F)
	colnames(g) = varnames
	check = c(check,varnames)
	g = cbind(frontpart,g,stringsAsFactors=F)
	save(g,file=paste0('chr',args[1],'_cut',k,'.RData')) # This is my output line
}

# Check we captured all variants in cuts
length(unique(check))==numvar
table(unique(check) %in% unique(colnames(geno[,7:ncol(geno)])))

