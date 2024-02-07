#' Sequence Alignment
#'
#' Performs sequence alignment on two sequences by a user specified alignment method. As of
#' now, only the Needleman-Wunsch algorithm is supported.
#'
#' @param seq1 First sequence to align.
#' @param seq2 Second sequence to align.
#' @param d Gap penalty, should be negative.
#' @param mismatch Mismatch penalty, should be negative.
#' @param match Match score, should be positive.
#' @param method Name of alignment algorithm. Currently only supports "needleman".
#'
#' @return Object of class \code{alignment} representing the alignment result.
#' This object can be utilized with the \code{\link{plot.alignment}} function to visualize
#' the alignment matrix and the \code{\link{print.alignment}} function to display alignments
#' in the console.
#'
#' @references
#' For more details on the Needleman-Wunsch algorithm, see the \href{https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm}{wikipedia page}.
#'
#' @examples
#' seq1 <- "GCATGCG"
#' seq2 <- "GATTACA"
#' # Run the Needleman-Wunsch algorithm
#' align_sequences(seq1, seq2, d = -1, mismatch = -1, match = 1)
#'
#' @export
#'
align_sequences <- function(seq1, seq2, d, mismatch, match, method ="needleman") {
  test_input <- function(value, type, name) {
    if ((type == "numeric" && (!is.numeric(value) || length(value) != 1)) ||
        (type == "character" && (!is.character(value) || length(value) != 1))) {
      stop(paste(name, "must be a single", type))
    }
  }

  # Test input
  test_input(seq1, "character", "seq1")
  test_input(seq2, "character", "seq2")
  test_input(d, "numeric", "d")
  test_input(mismatch, "numeric", "mismatch")
  test_input(match, "numeric", "match")

  # Check if d or mismatch is larger than 0 and match is larger than 1
  if (d > 0) warning("d (gap penalty) should be less than or equal to 0.")
  if (mismatch > 0) warning("mismatch (penalty) should be less than or equal to 0.")
  if (match < 0) warning("match (score) should be less than or equal to 1.")

  if (method=="needleman") {
    align_data <- needleman(seq1, seq2, d, mismatch, match)
    # More algorithms can be added here... if(method=="NAME"){}...
    } else {
    stop("Invalid alignment algorithm, currently supported: 'needleman'")
    }

  class(align_data) <- "alignment"
  return(align_data)
}
