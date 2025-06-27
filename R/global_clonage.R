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

# Fonction pour extraire la couleur hexadécimale des features GenBank
extract_color <- function(feature_lines) {
  color_line <- grep("SerialCloner_Color=", feature_lines, value = TRUE)
  if (length(color_line) > 0) {
    # Extraire la couleur hex (format &hXXXXXX)
    color_match <- regmatches(color_line[1], regexpr("&h[0-9A-Fa-f]{6}", color_line[1]))
    if (length(color_match) > 0) {
      # Convertir &hXXXXXX en #XXXXXX
      hex_color <- gsub("&h", "#", color_match)
      return(hex_color)
    }
  }
  return("#000000") # couleur par défaut
}

# Parser features GenBank avec vraies couleurs
parse_features <- function(features_lines) {
  feats <- list()
  current_feat <- NULL
  current_feat_lines <- character()

  for (line in features_lines) {
    if (grepl("^\\s{5}\\S", line)) {
      # Nouvelle feature (commence par 5 espaces suivis d'un mot)
      if (!is.null(current_feat)) {
        # Traiter la feature précédente
        current_feat$color <- extract_color(current_feat_lines)
        feats[[length(feats) + 1]] <- current_feat
      }

      # Nouvelle feature
      pos_match <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
      current_feat <- list(
        pos_raw = if (length(pos_match) > 0) pos_match else NA,
        color = "#000000",
        name = ""
      )
      current_feat_lines <- character()
    }

    if (!is.null(current_feat)) {
      current_feat_lines <- c(current_feat_lines, line)

      # Extraire le label
      if (grepl("/label=", line)) {
        label_match <- regmatches(line, regexpr('/label=("[^"]*"|[^\\s]*)', line))
        if (length(label_match) > 0) {
          feature_name <- gsub('/label=', '', label_match)
          feature_name <- gsub('"', '', feature_name)
          current_feat$name <- feature_name
        }
      } else if (grepl("/gene=", line)) {
        # Fallback : essayer aussi le champ /gene=
        gene_match <- regmatches(line, regexpr('/gene=("[^"]*"|[^\\s]*)', line))
        if (length(gene_match) > 0 && current_feat$name == "") {
          gene_name <- gsub('/gene=', '', gene_match)
          gene_name <- gsub('"', '', gene_name)
          current_feat$name <- gene_name
        }
      }
    }
  }

  # Ajouter la dernière feature
  if (!is.null(current_feat)) {
    current_feat$color <- extract_color(current_feat_lines)
    feats[[length(feats) + 1]] <- current_feat
  }

  feats
}

# Construire vecteur couleurs pour la séquence
get_color_map <- function(features_lines, aln_start, seq_length) {
  feats <- parse_features(features_lines)
  color_map <- rep(NA, seq_length)
  for (feat in feats) {
    pos_str <- feat$pos_raw
    if (is.null(pos_str) || is.na(pos_str)) next
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

# Générer la légende des couleurs
generate_legend <- function(features_lines) {
  feats <- parse_features(features_lines)
  if (length(feats) == 0) return("")

  legend_items <- character()
  for (feat in feats) {
    pos_str <- feat$pos_raw
    if (is.na(pos_str)) next
    bounds <- as.numeric(unlist(strsplit(pos_str, "\\.\\.")))
    if (length(bounds) != 2) next
    name_display <- if (feat$name != "") feat$name else "Feature sans nom"
    legend_items <- c(legend_items,
                      paste0("<span style='color:", feat$color, "; font-weight: bold;'>■</span> ",
                             name_display, " (", bounds[1], "-", bounds[2], ")"))
  }

  paste0("<div id='legend-fixed' style='position: sticky; top: 0; z-index: 1000; margin-bottom: 15px; padding: 10px; background: #f0f0f0; border-radius: 5px; border: 1px solid #ddd;'>",
         "<b>Légende des couleurs :</b><br>",
         paste(legend_items, collapse = "<br>"),
         "</div>")
}

# Générer la règle de numérotation corrigée
generate_ruler <- function(start_pos, seq_length) {
  ruler_chars <- rep(" ", seq_length)

  # Placer les numéros tous les 10 nucléotides
  for (i in seq_len(seq_length)) {
    pos <- start_pos + i - 1

    # Marquer chaque position multiple de 10
    if (pos %% 10 == 0) {
      num_str <- as.character(pos)

      # Calculer où placer le numéro (centré sur la position)
      num_len <- nchar(num_str)
      center_offset <- floor(num_len / 2)
      start_idx <- i - center_offset

      # Placer chaque chiffre du numéro
      for (j in 1:num_len) {
        target_pos <- start_idx + j - 1
        if (target_pos >= 1 && target_pos <= seq_length) {
          ruler_chars[target_pos] <- substr(num_str, j, j)
        }
      }
    }
  }

  # Ajouter des marqueurs pour les positions multiples de 5 (mais pas 10)
  for (i in seq_len(seq_length)) {
    pos <- start_pos + i - 1
    if (pos %% 5 == 0 && pos %% 10 != 0 && ruler_chars[i] == " ") {
      ruler_chars[i] <- "|"
    }
  }

  paste(ruler_chars, collapse = "")
}

# Fonction pour créer ligne HTML caractère par caractère en monospace, avec couleur et tooltip
htmlify_line <- function(line, colors = NULL, start_pos = NULL) {
  chars <- strsplit(line, "")[[1]]
  spans <- character(length(chars))
  for (i in seq_along(chars)) {
    style <- "font-family: Courier New, monospace;"
    if (!is.null(colors) && !is.na(colors[i])) {
      style <- paste0(style, "color:", colors[i], "; font-weight: bold;")
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
