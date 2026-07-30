#Functions used in the Michaelis Menten fitting app

# Load file function
# adapted from H Wickham Mastering Shiny, section 9.1.3

load_file <- function(aName, aPath, aSheet) {
  ext <- tools::file_ext(aName)
  switch(ext,
         xlsx = readxl::read_excel(aPath, aSheet, .name_repair = "universal", col_types = "numeric"),
         csv = vroom::vroom(aPath, delim = ",", show_col_types = FALSE, .name_repair = "universal"),
         tsv = vroom::vroom(aPath, delim = "\t", .name_repair = "universal"),
         txt = vroom::vroom(aPath, show_col_types = FALSE, .name_repair = "universal"),
         validate("Invalid file. Please upload a data file")
  )
}
# note
# Define Style
my_style <- function(x) {
  x |>
    adjust_colors(colors_discrete_candy) |>
    add_data_points(shape = 23, size = 3, fill = "lightblue1", color="steelblue") |>
    add_mean_dot(color = "red") |> 
    adjust_font(family = "verdana", face = "bold") |>
    adjust_size(width = 150, height = 100) |> 
    theme_tidyplot(fontsize = 16)
}

# Set global options
tidyplots_options(my_style = my_style)

#Linear transform plots
linPlot <- function(aDF, xvar, yvar) {
  aDF |>
    tidyplot(x = {{ xvar }}, y = {{ yvar }}) |>
    adjust_title("Linear transform",face = "bold", fontsize = 18) |>
    adjust_x_axis_title("X", face = "bold") |>
    adjust_y_axis_title("Y", face = "bold") |>
    add_curve_fit(method = "lm", se = FALSE, color = "purple4", linewidth = 0.75)
}

#Nonlinear curve fits
mmPlot <- function(aDF, xvar, yvar, Km, Vmax) {
  aDF |>
    tidyplot(x = {{ xvar }}, y = {{ yvar }}) |>
    adjust_title("Nonlinear plot",face = "bold", fontsize = 18) |> 
    add_reference_lines(
      x = Km,
      y = Vmax,
      linetype = "dotdash", 
      linewidth = 0.5
    ) |>
    adjust_x_axis_title("[S]", face = "bold") |>
    adjust_y_axis_title("V", face = "bold") |>
    add_curve_fit(method = "nls", formula = y ~ SSmicmen(x, Vm, K), se = FALSE, colour = "purple4", linewidth = 0.75)
}