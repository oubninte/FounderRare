#'  Generate a random color that is not too dark
#'
#' @param N : Number of code color to generate
#'
#' @return Color code of N colors
#' @export
#'
#' @examples generate_random_colors(4)
#'
#
generate_random_colors <- function(N) {
  colors <- vector("character", N)  # Create a character vector to store the colors
  rgb_colors=vector("integer", N)

  for (i in 1:N) {
    hsl_color <- ""  # Default value for hsl_color
    rgb_color <- c(0, 0, 0)  # Initialize rgb_color with a dummy value


    while (sum(rgb_color) <  150 ) {
      rgb_color <- sample(0:225, 3)   # Generate random RGB values between 0 and 200
      hsl_color <- rgb(rgb_color[1], rgb_color[2], rgb_color[3], maxColorValue = 255)  # Convert RGB to hexadecimal color
      hsl_color <- rgb(rgb_color[1], rgb_color[2], rgb_color[3], maxColorValue = 255)  # Convert RGB to hexadecimal color
      rgb_colors[i]=sum(rgb_color)
      if (i>=2 && (rgb_colors[i]-rgb_colors[i-1])<30) {
        rgb_color=rep(50,3)
      }
    }

    colors[i] <- hsl_color  # Store the generated color in the vector
  }

  return(colors)
}
