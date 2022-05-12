#!/usr/bin/env Rscript

script <- sub('--file=','',grep('file',commandArgs(),value=T))
if(length(script) > 0) {
	line <- paste0(rep("=",max(10,nchar(script)+1)), collapse="")
	text <- paste0("\n",line,"\nRscript:\n ",script," ",paste(commandArgs(T), collapse=" "),"\n",line,"\n")
	writeLines(text)
} 

