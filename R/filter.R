#'Filter.
#'
#' Filter low expressing counts by mean logarithm counts-per-million.
#'
#' @param data list consisting of count data and sample table.
#' @param threshold minimum mean-log2-CPM to keep the data.
#'
#' @import edgeR
#'
#' @return input list containing filtered data.
#' @export

filter <- function(data, threshold = 1){

#Extracts count data from input list:
count_data <- data[[1]]

#Filters count data:
mean_log2_cpm <- rowMeans(log2(cpm(count_data)+1))
count_data <- count_data[mean_log2_cpm > threshold,]

#Returns filtered data to input list:
data[[1]] <- count_data

return(data)
}
