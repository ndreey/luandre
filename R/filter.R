#'Filter.
#'
#' Filter low expressing counts by mean logarithm counts-per-million.
#'
#' @param data list consisting of the count data and sample table.
#' @param threshold minimum mean-log2-CPM to keep the data.
#'
#' @import edgeR
#'
#' @return input list containing filtered data.
#' @export

filter <- function(data, threshold = 1){

count_data <- data[[1]]

#Filters count data:
meanlog2cpm <- rowMeans(log2(cpm(count_data)+1))
count_data <- count_data[meanlog2cpm > threshold,]

data[[1]] <- count_data

return(data)
}
