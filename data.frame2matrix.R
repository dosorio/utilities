#' Converts three columns of a data.frame into a matrix -- e.g. to plot 
#' the data via image() later on. Two of the columns form the row and
#' col dimensions of the matrix. The third column provides values for
#' the matrix.
#' 
#' @param data data.frame: input data
#' @param rowtitle string: row-dimension; name of the column in data, which distinct values should be used as row names in the output matrix
#' @param coltitle string: col-dimension; name of the column in data, which distinct values should be used as column names in the output matrix
#' @param datatitle string: name of the column in data, which values should be filled into the output matrix
#' @param rowdecreasing logical: should the row names be in ascending (FALSE) or in descending (TRUE) order?
#' @param coldecreasing logical: should the col names be in ascending (FALSE) or in descending (TRUE) order?
#' @param default_value numeric: default value of matrix entries if no value exists in data.frame for the entries
#' @return matrix: matrix containing values of data[[datatitle]] with rownames data[[rowtitle]] and colnames data[coltitle]
#' @author Daniel Neumann
#' @date 2017-08-29
data.frame2matrix = function(data, rowtitle, coltitle, datatitle, 
                             rowdecreasing = FALSE, coldecreasing = FALSE,
                             default_value = NA) {

  # check, whether titles exist as columns names in the data.frame data
  if ( (!(rowtitle%in%names(data))) 
       || (!(coltitle%in%names(data))) 
       || (!(datatitle%in%names(data))) ) {
    stop('data.frame2matrix: bad row-, col-, or datatitle.')
  }

  # get number of rows in data
  ndata = dim(data)[1]

  # extract rownames and colnames for the matrix from the data.frame
  rownames = sort(unique(data[[rowtitle]]), decreasing = rowdecreasing)
  nrows = length(rownames)
  colnames = sort(unique(data[[coltitle]]), decreasing = coldecreasing)
  ncols = length(colnames)

  # initialize the matrix
  out_matrix = matrix(NA, 
                      nrow = nrows, ncol = ncols,
                      dimnames=list(rownames, colnames))

  # iterate rows of data
  for (i1 in 1:ndata) {
    # get matrix-row and matrix-column indices for the current data-row
    iR = which(rownames==data[[rowtitle]][i1])
    iC = which(colnames==data[[coltitle]][i1])

    # throw an error if the matrix entry (iR,iC) is already filled.
    if (!is.na(out_matrix[iR, iC])) stop('data.frame2matrix: double entry in data.frame')
    out_matrix[iR, iC] = data[[datatitle]][i1]
  }

  # set empty matrix entries to the default value
  out_matrix[is.na(out_matrix)] = default_value

  # return matrix
  return(out_matrix)

}
