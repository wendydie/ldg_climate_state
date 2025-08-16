#' Check whether hemispheres have adjacent palaeolatitudinal bins within 0–60°
#'
#' This function assesses whether each hemisphere (Northern or Southern)
#' in each time slice (`stage`) contains at least two adjacent palaeolatitudinal
#' bins within the tropical and temperate zone (0–60°). Hemispheres failing this
#' criterion are labeled as "bad"; others are labeled "good".
#'
#' @param rich_df A data frame containing at least `bin_midpoint`, `hemisphere`, and `bin` columns.
#' @param lat_bins A data frame defining the latitude bins, with `bin`, `min`, and `max` columns.
#'
#' @return A modified `rich_df` with an additional column `label` ("good" or "bad")
#'         indicating whether the hemisphere has sufficient spatial coverage for analysis.
#'
#' @examples
#' rich_df <- has_adjacent_bins(rich_df, lat_bins)
#'
#' @export
has_adjacent_bins <- function(rich_df, lat_bins) {
  # Define target bins for each hemisphere
  north_bins <- lat_bins$bin[lat_bins$min >= 0 & lat_bins$max <= 60]
  south_bins <- lat_bins$bin[lat_bins$min >= -60 & lat_bins$max <= 0]
  
  # Identify whether each hemisphere at each time slice has adjacent bins in 0–60°
  label_df <- rich_df %>%
    distinct(bin_midpoint, hemisphere, bin) %>%
    group_by(bin_midpoint, hemisphere) %>%
    summarise(
      label = {
        h <- unique(hemisphere)
        bins <- sort(unique(bin))
        target <- if (h == "Northern") north_bins
        else if (h == "Southern") south_bins
        else integer(0)
        b <- intersect(bins, target)
        if (length(b) >= 3 && any(diff(b) == 1)) "good" else "bad"
      },
      .groups = "drop"
    )
  
  
  # Merge label back into original dataframe
  left_join(rich_df, label_df, by = c("bin_midpoint", "hemisphere"))
}

