library("foreign")

files <- list.files(path="path/to/dir", pattern="*.txt", full.names=T, recursive=FALSE)
lapply(files, function(x) {
    	t <- read.table(x, header=T) # load file
    	# apply function
	mydata=read.csv(files)
	write.arff(x =mydata ,file= ".arff")
    	out <- function(t)
    	# write to file
    	write.table(out, "path/to/output", sep="\t", quote=F, row.names=F, col.names=T)
})
