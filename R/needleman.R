#' Initialize value matrix
#'
#' Creates a matrix with column names as \code{seq1} and
#' row names as \code{seq2}.
#'
#' @param seq1 First sequence to align.
#' @param seq2 Second sequence to align.
#' @param fill Value to fill the matrix with.
#'
#' @return A matrix with dimensions (\code{length(seq2)} x \code{length(seq1)}).
#' Filled with \code{fill} value.
#' @noRd
#'
init_mat <- function(seq1, seq2, fill){
  mat <- matrix(fill, ncol = length(seq1) + 1, nrow = length(seq2) + 1)
  colnames(mat) <- c("GAP", seq1)
  rownames(mat) <- c("GAP", seq2)
  return(mat)
}

#' Initialize direction list
#'
#' Creates a list of matrices. Can be used to keep track of the existence
#' of an arrow in a given location and direction.
#'
#' @param seq1 First sequence to align.
#' @param seq2 Second sequence to align.
#' @param fill Value to fill matrices with.
#'
#' @return A list of matrices, containing \code{up}, \code{left}, and \code{diagonal}
#' each matrix is filled with the value \code{fill}.
#' @noRd
#'
init_direction_lst <- function(seq1, seq2, fill=FALSE){
  arrow_mat <- init_mat(seq1,seq2, fill)
  dir_lst <- list("up" = arrow_mat, "left" = arrow_mat, "diagonal" = arrow_mat)
  return(dir_lst)
}

#' Define penalty
#'
#' @param d Gap penalty.
#' @param len_seq Sequence length.
#' @return Linear gap penalty.
#' @noRd
linear_penalty <- function(d, len_seq) {
  return(seq(0, len_seq-1) * d)
}

#' Get value in matrix
#'
#' @param mat value matrix
#' @param i row index
#' @param j column index
#' @param d Gap penalty
#' @param mismatch mismatch penalty
#' @param match score for a match
#' @param arrow_lst list of boolean arrow matrices
#'
#' @return list with the highest value and the updated arrow_lst
#' @references Following the Needleman-Wunsch algorithm (\href{https://en.wikipedia.org/wiki/Depth-first_search}{wikipedia page})
#' @noRd
#'
get_val <- function(mat, i, j, d, mismatch, match, arrow_lst){
  # TRUE for match, FALSE otherwise
  match_flag <- (colnames(mat)[j] == rownames(mat)[i])

  # Up
  up <- mat[i-1, j] + d
  # Left
  left <- mat[i, j-1] + d
  # Diagonal
  if (match_flag){
    dia <- mat[i-1, j-1] + match
  }  else {
    dia <- mat[i-1, j-1] + mismatch
  }

  max_val <- max(up, left, dia)
  ind_max <- which(c(up, left, dia)==max_val)

  for (index in ind_max){
    arrow_lst[[index]][i, j] <- TRUE
  }

  return(list(max_val, arrow_lst))
}

#' Compute optimal paths for Needleman-Wunsch algorithm
#'
#' Finds the optimal alignment by running Depth-First Search.
#'
#' @param arrow_lst List of matrices indicating the existence of arrows in different directions.
#'
#' @return List of directions to take from the lower right corner to the upper left corner.
#'
#' @details
#' The function uses Depth-First Search to find all optimal paths from the lower right corner
#' to the upper left corner in the alignment matrix.
#'
#' @noRd
#'
optimal_paths <- function(arrow_lst) {
  # Get starting location
  x <- ncol(arrow_lst$up)
  y <- nrow(arrow_lst$up)

  # Find neighboring cells
  neighbours <- function(arrow_lst, x, y) {
    result <- character(0)
    if (arrow_lst$diagonal[y, x]) {
      result <- c(result, "diagonal")
    }
    if (arrow_lst$left[y, x]) {
      result <- c(result, "left")
    }
    if (arrow_lst$up[y, x]) {
      result <- c(result, "up")
    }
    return(result)
  }

  # Depth-First Search to find all optimal paths
  dfs <- function(x, y, path, res_list) {
    if (x == 1 && y == 1) {
      res_list <- c(res_list, list(rev(path)))
      return(res_list)
    }

    nei <- neighbours(arrow_lst, x, y)
    for (neighbor in nei) {
      new_x <- x
      new_y <- y

      if (neighbor == "up") {
        new_y <- new_y - 1
      } else if (neighbor == "left") {
        new_x <- new_x - 1
      } else if (neighbor == "diagonal") {
        new_x <- new_x - 1
        new_y <- new_y - 1
      }

      res_list <- dfs(new_x, new_y, c(path, neighbor), res_list)
    }
    return(res_list)
  }

  res_list <- list()
  res_list <- dfs(x, y, character(0), res_list)

  # change order of result
  res_list <- lapply(res_list, function(path) rev(path))
  return(res_list)
}


#' Get alignment from paths for Needleman-Wunsch algorithm
#'
#' @param paths List returned by \code{\link{optimal_paths}}
#' @param seq1 First sequence to align.
#' @param seq2 Second sequence to align.
#'
#' @return A list containing:
#'   \item{path_loc}{Optimal path cells.}
#'   \item{alignments}{Aligned sequences.}
#'   \item{arrow_color}{Color of the arrows in the graph.}
#'
#' @details
#' The function processes the optimal paths obtained from the Needleman-Wunsch algorithm
#' and generates aligned sequences and additional information for visualization purposes.
#'
#' @examples
#' seq1 <- "GCATGCG"
#' seq2 <- "GATTACA"
#' paths <- optimal_paths(arrow_lst)
#' alignment_info <- get_alignment(paths, seq1, seq2)
#'
#' @noRd
#'
get_alignment <- function(paths, seq1, seq2){
  seq1_len <- length(seq1)
  seq2_len <- length(seq2)

  optimal_path_loc <- init_mat(seq1, seq2, fill=FALSE)
  optimal_path_loc[1,1] <- TRUE
  alignments <- list()

  arrow_col <- init_direction_lst(seq1, seq2, "gray")
  if (length(paths)>0) {
    for (i in 1:length(paths)) {
      path <- paths[[i]]
      seq1_alignment <- c()
      seq2_alignment <- c()
      seq1_rev <- rev(seq1)
      seq2_rev <- rev(seq2)
      x <- seq1_len + 1
      y <- seq2_len + 1

      for (j in 1:length(path)) {
        step <- path[j]
        optimal_path_loc[y,x] <- TRUE
        current_x_letter <- seq1[j]
        current_y_letter <- seq2[j]

        if (step=="diagonal") {
          seq1_letter <- seq1_rev[1]
          seq1_alignment[j] <- seq1_letter
          seq1_rev <- seq1_rev[-1]

          seq2_letter <- seq2_rev[1]
          seq2_alignment[j] <- seq2_letter
          seq2_rev <- seq2_rev[-1]

          if (seq1_letter==seq2_letter) {
            arrow_col$diagonal[y, x] <- "blue"
          } else {
            arrow_col$diagonal[y, x] <- "red"
          }

          x <- x-1
          y <- y-1
        }
        if (step=="up") {
          seq1_alignment[j] <- "-"
          seq2_alignment[j] <- seq2_rev[1]
          seq2_rev <- seq2_rev[-1]

          arrow_col$up[y, x] <- "#52e305"

          y <- y-1
        }
        if (step=="left") {
          seq1_alignment[j] <- seq1_rev[1]
          seq2_alignment[j] <- "-"
          seq1_rev <- seq1_rev[-1]

          arrow_col$left[y, x] <- "#52e305"

          x <- x-1
        }
      }
      alignments[[i]] <- list("Seq1"=paste(rev(seq1_alignment), collapse=""),
                              "Seq2"=paste(rev(seq2_alignment), collapse=""))
    }
  }
  return(list("path_loc"=optimal_path_loc,
              "alignments"=alignments,
              "arrow_color"=arrow_col))
}

#' Needleman-Wunsch Algorithm
#'
#' Performs global sequence alignment using the Needleman-Wunsch algorithm.
#' The function aligns two input sequences, allowing customization of gap
#' penalties, mismatch penalties, and match scores.
#'
#' @param seq1 First sequence to align.
#' @param seq2 Second sequence to align.
#' @param d Gap penalty, should be negative.
#' @param mismatch Mismatch penalty, should be negative.
#' @param match Match score, should be positive.
#'
#' @return A list containing:
#'     \item{arrows}{A list of matrices describing the arrows in the alignment graph.}
#'     \item{value_mat}{The alignment matrix.}
#'     \item{path_loc}{Cells part of the optimal alignment.}
#'     \item{alignments}{Aligned sequences.}
#'     \item{arrow_color}{Colors of the arrows.}
#'     \item{metadata}{A list of metadata from the function call.}
#'
#' @details
#' The Needleman-Wunsch algorithm is a dynamic programming algorithm used for global sequence alignment.
#' It is widely employed in bioinformatics to find the optimal alignment between two sequences,
#' allowing researchers to compare biological sequences such as DNA, RNA, or protein sequences.
#'
#' @examples
#' seq1 <- "GCATGCG"
#' seq2 <- "GATTACA"
#' # Run the Needleman-Wunsch algorithm
#' needleman(seq1, seq2, d = -1, mismatch = -1, match = 1)
#'
#' @references
#' For more details, see the \href{https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm}{wikipedia page}
#' for the Needleman-Wunsch algorithm.
#' @noRd
#'
needleman <- function(seq1, seq2, d, mismatch, match) {
  # Split sequences
  seq1 <- strsplit(seq1, "")[[1]]
  seq2 <- strsplit(seq2, "")[[1]]

  # Initialize matrix
  mat <- init_mat(seq1, seq2, fill=0)
  mat[, 1] <- linear_penalty(d, nrow(mat))
  mat[1, ] <- linear_penalty(d, ncol(mat))

  # Initialize arrow_lst
  arrow_lst <- init_direction_lst(seq1, seq2)
  arrow_lst$up[-1,1] <- TRUE
  arrow_lst$left[1,-1] <- TRUE

  # Calculate values in matrix
  for (i in 2:nrow(mat)) {
    for (j in 2:ncol(mat)) {
      res_list <- get_val(mat, i, j, d, mismatch, match, arrow_lst)
      mat[i, j] <- res_list[[1]]
      arrow_lst <- res_list[[2]]
    }
  }

  # Find optimal alignment path
  opt_paths <- optimal_paths(arrow_lst)

  # List of alignment
  align_lst <- get_alignment(opt_paths, seq1, seq2)

  # Collect metadata
  metadata <- list("seq1" = seq1, "seq2" = seq2, "d" = d,
                   "mismatch" = mismatch, "match" = match)

  res <- list("arrows" = arrow_lst,
              "value_mat" = mat,
              "path_loc" = align_lst$path_loc,
              "alignments" = align_lst$alignments,
              "arrow_col" = align_lst$arrow_col,
              "metadata" = metadata)
  return(res)
}

