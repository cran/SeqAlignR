#' Print Alignments
#'
#' Prints the alignments between \code{seq1} and \code{seq2} with the highest score.
#'
#' @param x Object of class \code{alignment} (see \code{\link{align_sequences}}).
#' @param ... Additional parameters to be passed to the \code{cat()} function, displaying the alignment.
#'
#' @examples
#' seq1 <- "GCATGCG"
#' seq2 <- "GATTACA"
#' # Run the Needleman-Wunsch algorithm
#' alignment1 <- align_sequences(seq1, seq2, d = -1, mismatch = -1, match = 1)
#' # Print the alignments
#' print(alignment1)
#'
#' @details
#' The printed message includes the alignment score.
#' This function may display multiple alignments, as alignments with the same score are possible.
#'
#' @return Console print of alignments.
#' @export
#'
print.alignment <- function(x, ...){
  alignments <- x$alignments
  mat <- x$value_mat
  md <- x$metadata
  x <- length(md$seq1) + 1
  y <- length(md$seq2) + 1
  cat("Alignments with a max score of ", mat[y,x],"\n", sep="")
  for (align in alignments) {
    s1 <- strsplit(align$Seq1, "")[[1]]
    s2 <- strsplit(align$Seq2, "")[[1]]
    space <- c()
    for(i in 1:length(s1)){
      if (s1[i]==s2[i]) {
        space <- c(space, "|")
      } else {
        space <- c(space, " ")
      }
    }
    cat(align$Seq1,"\n", ...)
    cat(paste(space, collapse=""), "\n", ...)
    cat(align$Seq2, ...)
    cat("\n\n")
  }
}
