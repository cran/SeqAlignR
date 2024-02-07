#' Plot Alignment Matrix
#'
#' Produces a plot displaying the alignment matrix of \code{seq1} and \code{seq2}.
#'
#' @param x Object of class \code{alignment} (see \code{\link{align_sequences}}).
#' @param ... Additional parameters to be passed to the \code{plot()} function.
#'
#' @examples
#' seq1 <- "GCATGCG"
#' seq2 <- "GATTACA"
#' # Run the Needleman-Wunsch algorithm
#' alignment1 <- align_sequences(seq1, seq2, d = -1, mismatch = -1, match = 1)
#' # Plot the matrix
#' plot(alignment1)
#'
#' @details
#' The first sequence (\code{seq1}) is represented by the columns and the second sequence (\code{seq2}) is represented
#' by the rows. The first column and first row are left bank, meaning a gap. Each cell in the matrix displays the score.
#' The subtitle states the \code{match}, \code{mismatch} and gap penalty \code{d} used in the algorithm.
#' A mismatch is shown by the red arrows, a match by the blue arrows, and a gap by the green arrows.
#' The alignment(s) with the highest score are highlighted with thick gray borders.
#'
#' @references
#' The implementation is inspired by the visualization \href{https://gist.github.com/slowkow/508393}{(code)} by Kamil Slowikowski (\href{https://orcid.org/0000-0002-2843-6370}{ORCID}).
#'
#' @return Plot of the alignment matrix.
#' @import graphics
#' @import plot.matrix
#' @export
#'
plot.alignment <- function(x, ...) {
  mat <- x$value_mat
  arrow_lst <- x$arrows
  metadata <- x$metadata
  path_loc <- x$path_loc
  arrow_col <- x$arrow_col

  oldpar <- par(mar = c(.8, .8, 4.8, .8), no.readonly = TRUE)
  on.exit(par(oldpar))
  colum_mat <- matrix(rep("", nrow(mat)), ncol=1)
  row_mat <- matrix(rep("", ncol(mat)+1), nrow=1)
  mat <- cbind(colum_mat, mat)
  mat <- rbind(row_mat, mat)
  plot(mat, fmt.cell = '%s', key = NULL, col = "white",
       ylab = "", xlab = "",
       axis.col = NULL, axis.row = NULL,
       main = paste0(""), ...)

  # Add bold labels for the first row
  text(1:ncol(mat), nrow(mat), labels = c("", "", colnames(mat)[3:ncol(mat)]), font = 2)

  # Add bold labels for the first column
  text(1, nrow(mat):1, labels = c("", "", rownames(mat)[3:nrow(mat)]), font = 2)

  mtext("Needleman-Wunsch", side = 3, line = 2.8, cex = 1.5, font = 1.9)

  mtext(paste("match =", metadata$match),
        side = 3, line = 0.8, adj=0.2,
        col = "blue", cex = 1.2)

  mtext(paste("mismatch =", metadata$mismatch),
        side = 3, line = 0.8, adj=0.5,
        col = "red", cex = 1.2)

  mtext(paste("gap =", metadata$d),
        side = 3, line = 0.8, adj=0.8,
        col = "#52e305", cex = 1.2)

  # change to gray/black border - area with/without arrows
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      i <- nrow(path_loc) - i + 2
      rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, lwd = 1,
           border = "#c1c1c1" # lightgray
      )
    }
  }

  # Add borders to mark solution(s)
  for (i in 1:nrow(path_loc)) {
    for (j in 1:ncol(path_loc)) {
      if (path_loc[i, j]) {
        rm <- nrow(path_loc) - (i + 1) + 2
        rect(j + 0.5, rm - 0.5, j + 1.5, rm + 0.5, lwd = 3.7,
             border = "#ababab",
        )
      }
    }
  }

  diff <- 0.34 # To adjust size of arrows. Smaller value -> longer arrows
  arrow_th <- 2.8 # Arrow thickness
  n_row <- nrow(arrow_lst$up)
  n_col <- ncol(arrow_lst$up)
  # Add diagonal arrows
  for (i in 2:nrow(arrow_lst$diagonal)) {
    for (j in 2:ncol(arrow_lst$diagonal)) {
      if (arrow_lst$diagonal[i, j]) {
        x1 <- j + diff
        y1 <- n_row - i - diff + 2
        x0 <- j - diff +1
        y0 <- n_row - i + diff + 1
        arrows(x0, y0, x1, y1, length=0.1,
               col=arrow_col$diagonal[i,j], lwd=arrow_th)
      }
    }
  }

  # Add left arrows
  for (i in 1:nrow(arrow_lst$left)) {
    for (j in 1:ncol(arrow_lst$left)) {
      if (arrow_lst$left[i, j]) {
        x1 <- j + diff
        y1 <- n_row - i + 1
        x0 <- j - diff + 1
        y0 <- n_row - i + 1
        arrows(x0, y0, x1, y1, length=0.1,
               col=arrow_col$left[i,j], lwd=arrow_th)
      }
    }
  }

  # Add up arrows
  for (i in 1:nrow(arrow_lst$up)) {
    for (j in 1:ncol(arrow_lst$up)) {
      if (arrow_lst$up[i, j]) {
        y1 <- n_row - i + 2 - diff
        x1 <- j +1
        y0 <- n_row - i + 1 + diff
        x0 <- j +1
        arrows(x0, y0, x1, y1, length=0.1,
               col=arrow_col$up[i,j], lwd=arrow_th)
      }
    }
  }
}

