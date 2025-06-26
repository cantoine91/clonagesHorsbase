library(Biostrings)

# Nettoyage séquence
clean_seq <- function(x) {
  x <- gsub("[\r\n ]", "", x)
  toupper(x)
}

# Annotation mutations (barre verticale)
annotate_mutations <- function(pattern_seq, subject_seq) {
  pattern_seq <- toupper(pattern_seq)
  subject_seq <- toupper(subject_seq)
  if (nchar(pattern_seq) != nchar(subject_seq)) {
    return(paste(rep(" ", nchar(pattern_seq)), collapse = ""))
  }
  pattern_chars <- strsplit(pattern_seq, "")[[1]]
  subject_chars <- strsplit(subject_seq, "")[[1]]
  annotation <- character(length(pattern_chars))
  for (i in seq_along(pattern_chars)) {
    annotation[i] <- ifelse(pattern_chars[i] != subject_chars[i], "|", " ")
  }
  paste(annotation, collapse = "")
}

# Parser features GenBank basique, extrait positions et couleur (exemple fixe)
parse_features <- function(features_lines) {
  feats <- list()
  for (line in features_lines) {
    pos_match <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
    if (length(pos_match) > 0) {
      feats[[length(feats) + 1]] <- list(pos_raw = pos_match, color = sample(c("red","blue","green","orange","purple"), 1))
    }
  }
  feats
}

# Construire vecteur couleurs pour la séquence
get_color_map <- function(features_lines, aln_start, seq_length) {
  feats <- parse_features(features_lines)
  color_map <- rep(NA, seq_length)
  for (feat in feats) {
    pos_str <- feat$pos_raw
    if (is.null(pos_str)) next
    bounds <- as.numeric(unlist(strsplit(pos_str, "\\.\\.")))
    if (length(bounds) != 2) next
    start <- bounds[1]
    end <- bounds[2]
    for (i in seq_len(seq_length)) {
      ref_pos <- aln_start + i - 1
      if (ref_pos >= start && ref_pos <= end) {
        color_map[i] <- feat$color
      }
    }
  }
  color_map
}

# Fonction pour créer ligne HTML caractère par caractère en monospace, avec couleur et tooltip
htmlify_line <- function(line, colors = NULL, start_pos = NULL) {
  chars <- strsplit(line, "")[[1]]
  spans <- character(length(chars))
  for (i in seq_along(chars)) {
    style <- "font-family: Courier New, monospace;"
    if (!is.null(colors) && !is.na(colors[i])) {
      style <- paste0(style, "color:", colors[i], ";")
    }
    tooltip <- if (!is.null(start_pos)) {
      paste0("Pos ", start_pos + i - 1)
    } else {
      ""
    }
    spans[i] <- paste0("<span title='", tooltip, "' style='", style, "'>", chars[i], "</span>")
  }
  paste0(spans, collapse = "")
}
