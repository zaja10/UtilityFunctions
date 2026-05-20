#' @title Match Germplasm Names with Genomic Matrix IDs
#' @description Evaluates whether a primary germplasm name or its synonyms exist in a set of valid genomic IDs.
#' @param name Primary germplasm name (string or factor)
#' @param syn_string Comma-separated string of synonyms
#' @param valid_ids Character vector of valid IDs in the genomic relationship matrix
#' @return The matched ID as a character string, or NA_character_ if no match is found
#' @importFrom stringr str_trim str_split
#' @export
match_germplasm <- function(name, syn_string, valid_ids) {
  name <- as.character(name)
  syn_string <- as.character(syn_string)
  
  if (name %in% valid_ids) {
    return(name)
  }
  
  if (!is.na(syn_string) && syn_string != "") {
    syn_list <- stringr::str_trim(unlist(stringr::str_split(syn_string, ",")))
    matched_syns <- syn_list[syn_list %in% valid_ids]
    
    if (length(matched_syns) > 0) {
      return(matched_syns[1])
    }
  }
  
  return(NA_character_)
}
