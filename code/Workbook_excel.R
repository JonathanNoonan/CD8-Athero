library(openxlsx)

save_supp_table <- function(
    wb_path,
    table_number,
    data,
    title,
    description,
    title_row = 1,
    desc_row  = 2,
    data_row  = 4,
    data_col  = 1,
    overwrite_workbook = TRUE
){
  # Construct sheet name automatically
  sheet_name <- paste0("Supplementary Table ", table_number)

  # Load workbook or create new
  if (file.exists(wb_path)) {
    wb <- loadWorkbook(wb_path)
  } else {
    wb <- createWorkbook()
  }

  # Add sheet if missing
  if (!(sheet_name %in% names(wb))) {
    addWorksheet(wb, sheet_name)
  }

  # Clear sheet first (overwrite everything while preserving tab order)
  # Method: write blank over the used region
  used <- getSheetDims(wb, sheet = sheet_name)
  if (!is.null(used)) {
    writeData(wb, sheet_name, x = data.frame(), startRow = 1, startCol = 1, colNames = FALSE)
  }

  # Write title
  writeData(
    wb, sheet = sheet_name,
    x = title,
    startRow = title_row, startCol = data_col
  )

  # Write description
  writeData(
    wb, sheet = sheet_name,
    x = description,
    startRow = desc_row, startCol = data_col
  )

  # Write data
  writeData(
    wb, sheet = sheet_name,
    x = data,
    startRow = data_row, startCol = data_col
  )

  # Save workbook
  saveWorkbook(wb, wb_path, overwrite = overwrite_workbook)

  message("✓ Written: ", sheet_name)
}
