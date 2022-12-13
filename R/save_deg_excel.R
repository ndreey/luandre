#' Save the Differentially Expressed Genes data.frame to excel.
#'
#' Creates a stylish excel file based on the DEG data.frame using openxlsx. The
#' file extension (.xslx) is automatically added to inputed filename.
#'
#' @param deg_df data.frame with differently expressed genes
#' @param filename string naming the excel file. (.xlsx is added automatically)
#'
#' @return filename.xlsx
#' @export
#'

save_deg_excel <- function(deg_df, filename) {

  # This can become an IF STATMENT to check if they added .xlsx or not...
  # The file extension is automatically added with placeholders.
  file <- "%s.xlsx"
  file <- sprintf(file, filename)

  # We want a column called SYMBOLS, using rowNAmes = TRUE wont give us colname.
  # Creates gene names column.
  deg_df$SYMBOL <- rownames(deg_df)

  # Reorder by column index so SYMBOL becomes first column
  deg_df <- deg_df[, c(6, 1:5)]

  # Rounds the logFC, logCPM and F column
  deg_df[,2:4] <- round(deg_df[, 2:4], 3)

  # We can style much more but this is a good start!
  # Creates a stylish heading in the excel sheet
  hs <- openxlsx::createStyle(textDecoration = "BOLD",
                              fontColour = "#FFFFFF",
                              fontSize = 12,
                              fontName = "Arial Narrow",
                              fgFill = "#4F80BD",
                              halign = "center",
                              valign = "center")

  # Create a custom style with bold formatting to add on the SYMBOL column
  bold_style <- openxlsx::createStyle(textDecoration = "BOLD",
                                      fontSize = 10,
                                      fontName = "Arial Narrow",
                                      halign = "center",
                                      valign = "center",
                                      border = "TopBottomLetRight",
                                      borderStyle = "thin")

  # Writes data to excel file
  openxlsx::write.xlsx(deg_df,
                       file,
                       sheetName = "luandre",
                       colNames = TRUE,
                       rowNames = FALSE,
                       tabColour = "grey",
                       borders = c("surrounding", "columns", "rows"),
                       headerStyle = hs)

  # Open the saved Excel file
  wb <- openxlsx::loadWorkbook(file)

  # Apply the bold_style to the cells in SYMBOL column of the first worksheet
  openxlsx::addStyle(wb,
                     sheet = 1,
                     cols = 1,
                     rows = 2:nrow(deg_df),
                     style = bold_style)

  # Save the updated Excel file
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)

}
