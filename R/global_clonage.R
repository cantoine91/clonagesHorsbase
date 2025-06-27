library(Biostrings)



# Nettoyage séquence
clean_seq <- function(x) {
  x <- gsub("[\r\n ]", "", x)
  toupper(x)
}

# Créer la séquence pattern complète avec gaps avant et après
create_full_pattern_sequence <- function(aligned_pattern, start_pos, total_length) {
  full_seq <- rep("-", total_length)

  # Placer la séquence alignée à partir de start_pos
  aligned_chars <- strsplit(aligned_pattern, "")[[1]]

  for (i in seq_along(aligned_chars)) {
    pos <- start_pos + i - 1
    if (pos <= total_length) {
      full_seq[pos] <- aligned_chars[i]
    }
  }

  return(paste(full_seq, collapse = ""))
}

# Créer l'annotation complète avec espaces avant et après
create_full_annotation_sequence <- function(aligned_annot, start_pos, total_length) {
  full_annot <- rep(" ", total_length)

  # Placer l'annotation alignée à partir de start_pos
  annot_chars <- strsplit(aligned_annot, "")[[1]]

  for (i in seq_along(annot_chars)) {
    pos <- start_pos + i - 1
    if (pos <= total_length) {
      full_annot[pos] <- annot_chars[i]
    }
  }

  return(paste(full_annot, collapse = ""))
}

# Annotation mutations - seulement | pour les correspondances exactes
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
    p_char <- pattern_chars[i]
    s_char <- subject_chars[i]

    # Correspondance exacte (pas de gaps)
    if (p_char == s_char && p_char != "-" && s_char != "-") {
      annotation[i] <- "|"
    }
    # Tout le reste (gaps, mismatches)
    else {
      annotation[i] <- " "
    }
  }

  return(paste(annotation, collapse = ""))
}

# Fonction AMÉLIORÉE pour extraire la couleur hexadécimale des features GenBank
extract_color <- function(feature_lines) {

  # Chercher la ligne contenant SerialCloner_Color (plus flexible)
  color_lines <- grep("SerialCloner_Color\\s*=", feature_lines, value = TRUE)

  if (length(color_lines) > 0) {
    # Extraire le code couleur &hXXXXXX (plus robuste)
    color_match <- regmatches(color_lines[1], regexpr("&h[0-9A-Fa-f]{6}", color_lines[1]))

    if (length(color_match) > 0) {
      # Convertir &hXXXXXX en #XXXXXX
      hex_color <- gsub("&h", "#", color_match)
      return(hex_color)
    } else {
      cat("Aucun code couleur valide trouvé dans:", color_lines[1], "\n")
    }
  } else {
    cat("Aucune ligne SerialCloner_Color trouvée\n")
  }
  return("#000000") # couleur par défaut
}

# Parser features GenBank AMÉLIORÉ avec debug et parsing robuste
parse_features <- function(features_lines) {

  feats <- list()
  current_feat <- NULL
  current_feat_lines <- character()

  for (line_num in seq_along(features_lines)) {
    line <- features_lines[line_num]

    # Détecter une nouvelle feature (commence par 5+ espaces puis le type)
    # Patterns possibles: misc_feature, CDS, rep_origin, polyA_signal, etc.
    if (grepl("^\\s{5,}[a-zA-Z_]+\\s+", line)) {

      # Traiter la feature précédente
      if (!is.null(current_feat)) {
        current_feat$color <- extract_color(current_feat_lines)
        feats[[length(feats) + 1]] <- current_feat
      }

      # Extraire le type de feature
      type_match <- regmatches(line, regexpr("^\\s+([a-zA-Z_]+)", line))
      feature_type <- if (length(type_match) > 0) {
        gsub("^\\s+", "", type_match)
      } else {
        "unknown"
      }

      # Extraire la position (plusieurs patterns possibles)
      pos_raw <- NA

      # Pattern 1: position simple (123..456)
      pos_match <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
      if (length(pos_match) > 0) {
        pos_raw <- pos_match
      } else {
        # Pattern 2: complement(123..456)
        comp_match <- regmatches(line, regexpr("complement\\(\\d+\\.\\.\\d+\\)", line))
        if (length(comp_match) > 0) {
          pos_match <- regmatches(comp_match, regexpr("\\d+\\.\\.\\d+", comp_match))
          if (length(pos_match) > 0) {
            pos_raw <- pos_match
          }
        }
      }

      current_feat <- list(
        type = feature_type,
        pos_raw = pos_raw,
        color = "#000000",
        name = ""
      )
      current_feat_lines <- character()
    }

    if (!is.null(current_feat)) {
      current_feat_lines <- c(current_feat_lines, line)

      # Extraire le label (plusieurs patterns possibles)
      if (grepl("/label\\s*=", line)) {
        # Pattern 1: /label="nom"
        if (grepl('/label\\s*=\\s*"[^"]*"', line)) {
          label_match <- regmatches(line, regexpr('/label\\s*=\\s*"[^"]*"', line))
          feature_name <- gsub('/label\\s*=\\s*"([^"]*)"', '\\1', label_match)
        }
        # Pattern 2: /label=nom (sans guillemets)
        else {
          label_match <- regmatches(line, regexpr('/label\\s*=\\s*[^\\s]+', line))
          feature_name <- gsub('/label\\s*=\\s*', '', label_match)
        }

        if (exists("feature_name") && nchar(feature_name) > 0) {
          current_feat$name <- feature_name
        }
      }

      # Extraire gene comme fallback
      else if (grepl("/gene\\s*=", line) && current_feat$name == "") {
        if (grepl('/gene\\s*=\\s*"[^"]*"', line)) {
          gene_match <- regmatches(line, regexpr('/gene\\s*=\\s*"[^"]*"', line))
          gene_name <- gsub('/gene\\s*=\\s*"([^"]*)"', '\\1', gene_match)
        } else {
          gene_match <- regmatches(line, regexpr('/gene\\s*=\\s*[^\\s]+', line))
          gene_name <- gsub('/gene\\s*=\\s*', '', gene_match)
        }

        if (exists("gene_name") && nchar(gene_name) > 0) {
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

  for (i in seq_along(feats)) {
    feat <- feats[[i]]
  }

  return(feats)
}

# Configuration des couleurs par défaut (fallback)
default_feature_colors <- list(
  "rabbit IgG Fc [3']" = "#FF8000",
  "insert from Sequence Window #2" = "#84A4C0",
  "ColE1 origin" = "#4682B4",
  "ZeoR" = "#8FBC8F",
  "SV40 late polyA" = "#FFC0CB",
  "IL2 ss" = "#8FBC8F"
)

# Fonction pour obtenir la couleur par nom (fallback)
get_color_by_name <- function(feature_name) {
  color <- default_feature_colors[[feature_name]]
  if (is.null(color)) {
    return("#888888") # Gris par défaut
  }
  return(color)
}

# Construire vecteur couleurs pour la séquence AMÉLIORÉ
get_color_map <- function(features_lines, aln_start, seq_length) {

  feats <- parse_features(features_lines)
  color_map <- rep(NA, seq_length)

  for (feat in feats) {
    pos_str <- feat$pos_raw
    if (is.null(pos_str) || is.na(pos_str)) {
      next
    }

    bounds <- as.numeric(unlist(strsplit(pos_str, "\\.\\.")))
    if (length(bounds) != 2) {
      next
    }

    start_feat <- bounds[1]
    end_feat <- bounds[2]

    # Utiliser la couleur extraite ou le fallback
    color_to_use <- feat$color
    if (color_to_use == "#000000") {
      color_to_use <- get_color_by_name(feat$name)
    }
    # Appliquer la couleur aux positions correspondantes
    for (i in seq_len(seq_length)) {
      ref_pos <- aln_start + i - 1
      if (ref_pos >= start_feat && ref_pos <= end_feat) {
        color_map[i] <- color_to_use
      }
    }
  }

  # Compter les positions colorées
  colored_positions <- sum(!is.na(color_map))
  return(color_map)
}

# Générer la légende des couleurs AMÉLIORÉE
generate_legend <- function(features_lines) {
  feats <- parse_features(features_lines)
  if (length(feats) == 0) {
    return("<div id='legend-fixed'><b>Aucune annotation de couleur trouvée</b></div>")
  }

  legend_items <- character()
  for (feat in feats) {
    pos_str <- feat$pos_raw
    if (is.na(pos_str)) next

    bounds <- as.numeric(unlist(strsplit(pos_str, "\\.\\.")))
    if (length(bounds) != 2) next

    name_display <- if (feat$name != "") feat$name else paste("Feature", feat$type)

    # Utiliser la couleur extraite ou le fallback
    color_to_use <- feat$color
    if (color_to_use == "#000000") {
      color_to_use <- get_color_by_name(feat$name)
    }

    legend_items <- c(legend_items,
                      paste0("<span style='color:", color_to_use, "; font-weight: bold;'>■</span> ",
                             name_display, " (", bounds[1], "-", bounds[2], ")"))
  }

  result <- paste0("<div id='legend-fixed'>",
                   "<b>Légende des couleurs :</b><br>",
                   paste(legend_items, collapse = "<br>"),
                   "</div>")

  return(result)
}

# Générer un alignement avec table HTML pour alignement parfait
generate_blast_style_alignment <- function(pattern_seq, subject_seq, annot_seq, start_pos, colors, filename) {
  line_length <- 60
  seq_length <- nchar(subject_seq)

  html_lines <- character()
  text_lines <- character()

  html_lines <- c(html_lines, paste0("<div class='alignment-block'>"))
  html_lines <- c(html_lines, paste0("<div class='alignment-title'>Alignement avec ", filename, "</div>"))

  text_lines <- c(text_lines, paste0("=== Alignement avec ", filename, " ==="))

  for (start in seq(1, seq_length, by = line_length)) {
    end <- min(start + line_length - 1, seq_length)

    # Extraire les segments
    seq1_segment <- substr(subject_seq, start, end)  # Référence
    seq2_segment <- substr(pattern_seq, start, end)  # Query
    annot_segment <- substr(annot_seq, start, end)

    # Positions
    seq1_start <- start
    seq1_end <- end

    seq2_bases_before <- nchar(gsub("-", "", substr(pattern_seq, 1, start - 1)))
    seq2_bases_in_segment <- nchar(gsub("-", "", seq2_segment))
    seq2_start <- seq2_bases_before + 1
    seq2_end <- seq2_bases_before + seq2_bases_in_segment

    if (seq2_bases_in_segment == 0) {
      seq2_end <- seq2_bases_before
    }

    # Couleurs pour la référence (segment spécifique)
    segment_colors <- if (!is.null(colors) && length(colors) >= end) {
      colors[start:end]
    } else {
      NULL
    }

    # Créer une table HTML pour alignement parfait
    table_html <- create_alignment_table(seq1_segment, seq2_segment, annot_segment,
                                         seq1_start, seq1_end, seq2_start, seq2_end,
                                         segment_colors)

    html_lines <- c(html_lines, table_html)

    # Version texte pour export
    text_lines <- c(text_lines, "")
    text_lines <- c(text_lines, sprintf("Seq_1 %6d %s %d", seq1_start, seq1_segment, seq1_end))
    text_lines <- c(text_lines, sprintf("       %s", annot_segment))
    text_lines <- c(text_lines, sprintf("Seq_2 %6d %s %d", seq2_start, seq2_segment, seq2_end))
  }

  html_lines <- c(html_lines, "</div>")

  return(list(
    html = paste(html_lines, collapse = "\n"),
    text = paste(text_lines, collapse = "\n")
  ))
}

# Créer une table HTML pour alignement parfait - VERSION AMÉLIORÉE
create_alignment_table <- function(seq1_seg, seq2_seg, annot_seg, seq1_start, seq1_end, seq2_start, seq2_end, segment_colors) {
  # Diviser en caractères individuels
  seq1_chars <- strsplit(seq1_seg, "")[[1]]
  seq2_chars <- strsplit(seq2_seg, "")[[1]]
  annot_chars <- strsplit(annot_seg, "")[[1]]

  # Créer les cellules colorées pour seq1
  seq1_cells <- character()
  for (i in seq_along(seq1_chars)) {
    style <- "font-family: 'Courier New', monospace; text-align: center; padding: 0; margin: 0; border: 0; width: 1ch; min-width: 1ch; max-width: 1ch;"

    # Appliquer la couleur si disponible
    if (!is.null(segment_colors) && i <= length(segment_colors) && !is.na(segment_colors[i])) {
      style <- paste0(style, " color: ", segment_colors[i], "; font-weight: bold;")
    }

    seq1_cells <- c(seq1_cells, paste0("<td style='", style, "'>", seq1_chars[i], "</td>"))
  }

  # Créer les cellules pour l'annotation
  annot_cells <- character()
  for (i in seq_along(annot_chars)) {
    style <- "font-family: 'Courier New', monospace; text-align: center; padding: 0; margin: 0; border: 0; width: 1ch; min-width: 1ch; max-width: 1ch;"
    annot_cells <- c(annot_cells, paste0("<td style='", style, "'>", annot_chars[i], "</td>"))
  }

  # Créer les cellules pour seq2
  seq2_cells <- character()
  for (i in seq_along(seq2_chars)) {
    style <- "font-family: 'Courier New', monospace; text-align: center; padding: 0; margin: 0; border: 0; width: 1ch; min-width: 1ch; max-width: 1ch;"
    seq2_cells <- c(seq2_cells, paste0("<td style='", style, "'>", seq2_chars[i], "</td>"))
  }

  # Construire la table avec espaces intégrés dans le texte
  table_style <- "border-collapse: collapse; margin: 10px 0; font-family: 'Courier New', monospace;"

  # Préfixes avec largeur fixe
  prefix_style <- "font-family: 'Courier New', monospace; padding: 0; margin: 0; text-align: left; width: 100px; min-width: 100px; max-width: 100px; white-space: nowrap;"
  suffix_style <- "font-family: 'Courier New', monospace; padding: 0; margin: 0; text-align: left; width: 60px; min-width: 60px; max-width: 60px;"

  table_html <- paste0(
    "<table style='", table_style, "'>",
    "<tr>",
    "<td style='", prefix_style, "'>Seq_1 ", sprintf("%6d", seq1_start), " </td>",
    paste(seq1_cells, collapse = ""),
    "<td style='", suffix_style, "'> ", seq1_end, "</td>",
    "</tr>",
    "<tr>",
    "<td style='", prefix_style, "'></td>",
    paste(annot_cells, collapse = ""),
    "<td style='", suffix_style, "'></td>",
    "</tr>",
    "<tr>",
    "<td style='", prefix_style, "'>Seq_2 ", sprintf("%6d", seq2_start), " </td>",
    paste(seq2_cells, collapse = ""),
    "<td style='", suffix_style, "'> ", seq2_end, "</td>",
    "</tr>",
    "</table>"
  )

  return(table_html)
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
