table2pdf <- function(tab, file) {
  texfile <- paste(file, ".tex", sep="")
  cat("\\documentclass{article}\n\\begin{document}\n", file=texfile)
  print(xtable(tab), include.rownames=FALSE, floating=FALSE, 
        file=texfile, append=TRUE)
  cat("\\end{document}\n", file=texfile, append=TRUE)
  system(paste("pdflatex", texfile))
}
