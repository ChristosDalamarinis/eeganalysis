#' Get BioSemi Channel Information
#'
#' Returns information about standard BioSemi channel layouts
#' (64-channel, 128-channel, 256-channel systems)
#'
#' @param n_channels Number of channels in the recording
#'
#' @return Data frame with columns: channel_number, name, position_type
#'
#' @details
#' BioSemi systems typically use standard 10-20/10-10 electrode positions.
#' This function provides information about standard configurations.
#'
#' @examples
#' \dontrun{
#'   biosemi_64_info <- get_biosemi_channel_info(64)
#' }
#'
#' @export
get_biosemi_channel_info <- function(n_channels) {
  
  # Standard 64-channel BioSemi layout (10-20 system)
  if (n_channels == 64) {
    channels_64 <- c(
      # Frontal
      "Fp1", "Fpz", "Fp2",
      # Anterior frontal
      "AF7", "AF3", "AFz", "AF4", "AF8",
      # Frontal
      "F7", "F5", "F3", "F1", "Fz", "F2", "F4", "F6", "F8",
      # Fronto-temporal
      "FT7", "FT8",
      # Fronto-central
      "FC5", "FC3", "FC1", "FCz", "FC2", "FC4", "FC6",
      # Central
      "T7", "C5", "C3", "C1", "Cz", "C2", "C4", "C6", "T8",
      # Centro-temporal
      "TP7", "TP8",
      # Centro-parietal
      "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6",
      # Parietal
      "P7", "P5", "P3", "P1", "Pz", "P2", "P4", "P6", "P8",
      # Parieto-temporal
      "PO7", "PO8",
      # Parieto-occipital
      "PO3", "POz", "PO4",
      # Occipital
      "O1", "Oz", "O2"
    )
    
    return(data.frame(
      channel_number = 1:64,
      name = channels_64,
      position_type = "10-20 system",
      stringsAsFactors = FALSE
    ))
  }
  
  # For other numbers of channels, return generic info
  return(data.frame(
    channel_number = 1:n_channels,
    name = paste0("Ch", seq_len(n_channels)),
    position_type = "Unknown",
    stringsAsFactors = FALSE
  ))
}
