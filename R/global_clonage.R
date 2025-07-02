library(Biostrings)

# ==================== FONCTIONS UTILITAIRES ====================

# Nettoyage de séquence
clean_sequence <- function(sequence_text) {
  sequence_text <- gsub("[\r\n ]", "", sequence_text)
  toupper(sequence_text)
}

# ==================== FONCTIONS SITES DE RESTRICTION ====================

# Base de données des enzymes de restriction communes
get_restriction_enzymes <- function() {
  return(list(
    "EcoRI" = "GAATTC",
    "BamHI" = "GGATCC",
    "HindIII" = "AAGCTT",
    "XhoI" = "CTCGAG",
    "SalI" = "GTCGAC",
    "NotI" = "GCGGCCGC",
    "XbaI" = "TCTAGA",
    "SpeI" = "ACTAGT",
    "KpnI" = "GGTACC",
    "SacI" = "GAGCTC",
    "BglII" = "AGATCT",
    "PstI" = "CTGCAG",
    "SmaI" = "CCCGGG",
    "ApaI" = "GGGCCC",
    "ClaI" = "ATCGAT",
    "NdeI" = "CATATG",
    "NcoI" = "CCATGG",
    "BspEI" = "TCCGGA",
    "AflII" = "CTTAAG",
    "MluI" = "ACGCGT",
    "AscI" = "GGCGCGCC",
    "PacI" = "TTAATTAA",
    "SbfI" = "CCTGCAGG",
    "AvrII" = "CCTAGG",
    "NheI" = "GCTAGC",
    "BssHII" = "GCGCGC",
    "PmeI" = "GTTTAAAC",
    "SrfI" = "GCCCGGGC",
    "HpaI" = "GTTAAC",
    "ScaI" = "AGTACT"
  ))
}

# Trouver les sites d'une enzyme
find_restriction_sites <- function(sequence, enzyme_sequence) {
  if (is.null(sequence) || is.null(enzyme_sequence) || enzyme_sequence == "") {
    return(integer(0))
  }

  if (class(sequence)[1] == "DNAString") {
    sequence <- as.character(sequence)
  }

  sequence <- toupper(sequence)
  enzyme_sequence <- toupper(enzyme_sequence)

  # Utiliser gregexpr pour trouver TOUTES les occurrences d'un coup
  matches <- gregexpr(enzyme_sequence, sequence, fixed = TRUE)[[1]]

  # Si aucune correspondance trouvée
  if (matches[1] == -1) {
    return(integer(0))
  }

  # Retourner toutes les positions
  return(as.integer(matches))
}

# Créer une carte des positions de sites de restriction
build_restriction_position_map <- function(sequence_length, restriction_sites, alignment_start = 1) {
  is_restriction <- rep(FALSE, sequence_length)

  if (!is.null(restriction_sites) && length(restriction_sites) > 0) {
    enzymes <- get_restriction_enzymes()

    for (enzyme_name in names(restriction_sites)) {
      sites <- restriction_sites[[enzyme_name]]
      if (length(sites) > 0) {
        site_length <- nchar(enzymes[[enzyme_name]])

        for (site_pos in sites) {
          for (i in 0:(site_length - 1)) {
            pos <- site_pos + i
            if (pos >= alignment_start && pos <= (alignment_start + sequence_length - 1)) {
              map_pos <- pos - alignment_start + 1
              if (map_pos > 0 && map_pos <= sequence_length) {
                is_restriction[map_pos] <- TRUE
              }
            }
          }
        }
      }
    }
  }

  return(is_restriction)
}

# Combiner les couleurs des features et des sites de restriction
combine_color_maps <- function(feature_colors, restriction_colors, sequence_length) {
  combined_colors <- rep(NA, sequence_length)

  for (i in 1:sequence_length) {
    if (!is.na(restriction_colors[i])) {
      # Priorité aux sites de restriction
      combined_colors[i] <- restriction_colors[i]
    } else if (!is.na(feature_colors[i])) {
      combined_colors[i] <- feature_colors[i]
    }
  }

  return(combined_colors)
}

# Générer la légende pour les sites de restriction
generate_restriction_legend_formatted <- function(restriction_sites_list) {
  if (length(restriction_sites_list) == 0) {
    return("")
  }

  colors <- c("#FF0000", "#0000FF", "#00FF00", "#FF8000", "#8000FF", "#00FFFF")
  legend_items <- character()

  site_index <- 1
  for (enzyme_name in names(restriction_sites_list)) {
    sites <- restriction_sites_list[[enzyme_name]]
    if (length(sites) > 0) {
      color <- colors[((site_index - 1) %% length(colors)) + 1]
      enzymes <- get_restriction_enzymes()
      enzyme_seq <- enzymes[[enzyme_name]]

      sites_text <- paste(sites, collapse = ", ")
      legend_items <- c(legend_items,
                        paste0("<div style='margin-bottom: 5px;'><span style='color:", color, "; font-weight: bold; font-size: 16px;'>■</span> ",
                               "<span style='font-size: 12px;'>", enzyme_name, " (", enzyme_seq, ") - Sites: ", sites_text, "</span></div>"))
      site_index <- site_index + 1
    }
  }

  if (length(legend_items) > 0) {
    return(paste0("<h4 style='margin-top: 15px; color: #b22222;'>Sites de restriction</h4>",
                  paste(legend_items, collapse = "")))
  }

  return("")
}

# ==================== FONCTIONS D'ALIGNEMENT ====================

# Créer la séquence pattern complète avec gaps
create_full_pattern_with_gaps <- function(aligned_pattern, start_position, total_length) {
  full_sequence <- rep("-", total_length)
  aligned_chars <- strsplit(aligned_pattern, "")[[1]]

  for (i in seq_along(aligned_chars)) {
    position <- start_position + i - 1
    if (position <= total_length) {
      full_sequence[position] <- aligned_chars[i]
    }
  }

  return(paste(full_sequence, collapse = ""))
}

# Créer l'annotation complète avec espaces
create_full_annotation_with_spaces <- function(aligned_annotation, start_position, total_length) {
  full_annotation <- rep(" ", total_length)
  annotation_chars <- strsplit(aligned_annotation, "")[[1]]

  for (i in seq_along(annotation_chars)) {
    position <- start_position + i - 1
    if (position <= total_length) {
      full_annotation[position] <- annotation_chars[i]
    }
  }

  return(paste(full_annotation, collapse = ""))
}

# Annoter les mutations - seulement | pour les correspondances exactes
annotate_sequence_mutations <- function(pattern_sequence, subject_sequence) {
  pattern_sequence <- toupper(pattern_sequence)
  subject_sequence <- toupper(subject_sequence)

  if (nchar(pattern_sequence) != nchar(subject_sequence)) {
    return(paste(rep(" ", nchar(pattern_sequence)), collapse = ""))
  }

  pattern_chars <- strsplit(pattern_sequence, "")[[1]]
  subject_chars <- strsplit(subject_sequence, "")[[1]]
  annotation <- character(length(pattern_chars))

  for (i in seq_along(pattern_chars)) {
    if (pattern_chars[i] == subject_chars[i] && pattern_chars[i] != "-" && subject_chars[i] != "-") {
      annotation[i] <- "|"
    } else {
      annotation[i] <- " "
    }
  }

  return(paste(annotation, collapse = ""))
}

# ==================== EXTRACTION DES COULEURS SERIALCLONER ====================

# Extraire la couleur hexadécimale des features GenBank
extract_serialcloner_color <- function(feature_lines) {
  color_lines <- grep("SerialCloner_Color\\s*=", feature_lines, value = TRUE)

  if (length(color_lines) > 0) {
    color_match <- regmatches(color_lines[1], regexpr("&h[0-9A-Fa-f]{6}", color_lines[1]))
    if (length(color_match) > 0) {
      return(gsub("&h", "#", color_match))
    }
  }

  return("#000000")
}

# ==================== PARSING DES FEATURES GENBANK ====================

# Parser les features GenBank de manière robuste
parse_genbank_features <- function(features_lines) {
  features_list <- list()
  current_feature <- NULL
  current_feature_lines <- character()

  for (i in seq_along(features_lines)) {
    line <- features_lines[i]
    is_new_feature_with_type <- grepl("^\\s{5,}[a-zA-Z_]+\\s+", line)
    is_position_only <- grepl("^\\s+\\d+\\.\\.\\d+", line)

    if (is_new_feature_with_type || is_position_only) {
      if (!is.null(current_feature)) {
        current_feature$color <- extract_serialcloner_color(current_feature_lines)
        features_list[[length(features_list) + 1]] <- current_feature
      }

      feature_type <- if (is_new_feature_with_type) extract_feature_type(line) else "misc_feature"
      position_raw <- extract_feature_position(line)

      current_feature <- list(
        type = feature_type,
        position_raw = position_raw,
        color = "#000000",
        name = ""
      )
      current_feature_lines <- character()

      # Pour les positions orphelines, regarder la ligne suivante pour le label
      if (is_position_only && i < length(features_lines)) {
        next_line <- features_lines[i + 1]
        if (grepl("/label\\s*=", next_line)) {
          feature_name <- extract_feature_name(next_line)
          if (!is.null(feature_name)) {
            current_feature$name <- feature_name
          }
        }
      }
    }

    if (!is.null(current_feature)) {
      current_feature_lines <- c(current_feature_lines, line)
      feature_name <- extract_feature_name(line)
      if (!is.null(feature_name) && current_feature$name == "") {
        current_feature$name <- feature_name
      }
    }
  }

  if (!is.null(current_feature)) {
    current_feature$color <- extract_serialcloner_color(current_feature_lines)
    features_list[[length(features_list) + 1]] <- current_feature
  }

  return(features_list)
}

# Extraire le type de feature
extract_feature_type <- function(line) {
  type_match <- regmatches(line, regexpr("^\\s+([a-zA-Z_]+)", line))
  if (length(type_match) > 0) {
    return(gsub("^\\s+", "", type_match))
  }
  return("unknown")
}

# Extraire la position de la feature
extract_feature_position <- function(line) {
  position_match <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
  if (length(position_match) > 0) {
    return(position_match)
  }

  complement_match <- regmatches(line, regexpr("complement\\(\\d+\\.\\.\\d+\\)", line))
  if (length(complement_match) > 0) {
    position_match <- regmatches(complement_match, regexpr("\\d+\\.\\.\\d+", complement_match))
    if (length(position_match) > 0) {
      return(position_match)
    }
  }

  return(NA)
}

# Extraire le nom de la feature (label ou gene)
extract_feature_name <- function(line) {
  if (grepl("/label\\s*=", line)) {
    if (grepl('/label\\s*=\\s*"[^"]*"', line)) {
      label_match <- regmatches(line, regexpr('/label\\s*=\\s*"[^"]*"', line))
      return(gsub('/label\\s*=\\s*"([^"]*)"', '\\1', label_match))
    } else {
      # Pour les labels sans guillemets, prendre tout ce qui suit =
      label_match <- regmatches(line, regexpr('/label\\s*=\\s*.+', line))
      if (length(label_match) > 0) {
        return(gsub('/label\\s*=\\s*', '', label_match))
      }
    }
  }

  if (grepl("/gene\\s*=", line)) {
    if (grepl('/gene\\s*=\\s*"[^"]*"', line)) {
      gene_match <- regmatches(line, regexpr('/gene\\s*=\\s*"[^"]*"', line))
      return(gsub('/gene\\s*=\\s*"([^"]*)"', '\\1', gene_match))
    } else {
      gene_match <- regmatches(line, regexpr('/gene\\s*=\\s*.+', line))
      if (length(gene_match) > 0) {
        return(gsub('/gene\\s*=\\s*', '', gene_match))
      }
    }
  }

  return(NULL)
}

# ==================== COULEURS PAR DÉFAUT ====================

# Configuration des couleurs par défaut pour les features connues
get_default_feature_colors <- function() {
  return(list(
    "rabbit IgG Fc [3']" = "#FF8000",
    "insert from Sequence Window #2" = "#84A4C0",
    "ColE1 origin" = "#4682B4",
    "ZeoR" = "#8FBC8F",
    "SV40 late polyA" = "#FFC0CB",
    "IL2 ss" = "#8FBC8F"
  ))
}

# Obtenir la couleur par nom de feature
get_color_by_feature_name <- function(feature_name) {
  default_colors <- get_default_feature_colors()
  color <- default_colors[[feature_name]]
  if (is.null(color)) {
    return("#888888")
  }
  return(color)
}

# ==================== MAPPING DES COULEURS ====================

# Construire le vecteur de couleurs pour la séquence
build_sequence_color_map <- function(features_lines, alignment_start, sequence_length) {
  features_list <- parse_genbank_features(features_lines)
  color_map <- rep(NA, sequence_length)

  for (feature in features_list) {
    position_string <- feature$position_raw
    if (is.null(position_string) || is.na(position_string)) {
      next
    }

    position_bounds <- as.numeric(unlist(strsplit(position_string, "\\.\\.")))
    if (length(position_bounds) != 2) {
      next
    }

    feature_start <- position_bounds[1]
    feature_end <- position_bounds[2]

    color_to_use <- feature$color
    if (color_to_use == "#000000") {
      color_to_use <- get_color_by_feature_name(feature$name)
    }

    for (i in seq_len(sequence_length)) {
      reference_position <- alignment_start + i - 1
      if (reference_position >= feature_start && reference_position <= feature_end) {
        color_map[i] <- color_to_use
      }
    }
  }

  return(color_map)
}

# ==================== GÉNÉRATION DE LA LÉGENDE ====================

# Générer la légende HTML des couleurs
generate_color_legend <- function(features_lines, restriction_sites_list = NULL) {
  features_list <- parse_genbank_features(features_lines)
  legend_items <- character()

  # Légende des features
  for (feature in features_list) {
    position_string <- feature$position_raw
    if (is.na(position_string)) next

    position_bounds <- as.numeric(unlist(strsplit(position_string, "\\.\\.")))
    if (length(position_bounds) != 2) next

    display_name <- if (feature$name != "") feature$name else paste("Feature", feature$type)

    color_to_use <- feature$color
    if (color_to_use == "#000000") {
      color_to_use <- get_color_by_feature_name(feature$name)
    }

    legend_items <- c(legend_items,
                      paste0("<div style='margin-bottom: 5px;'><span style='color:", color_to_use, "; font-weight: bold; font-size: 16px;'>■</span> ",
                             "<span style='font-size: 12px;'>", display_name, " (", position_bounds[1], "-", position_bounds[2], ")</span></div>"))
  }

  legend_content <- if (length(legend_items) > 0) {
    paste0("<h4 style='margin-top: 0; color: #b22222;'>Légende des couleurs</h4>",
           paste(legend_items, collapse = ""))
  } else {
    "<h4 style='margin-top: 0; color: #b22222;'>Aucune annotation trouvée</h4>"
  }

  # Ajouter la légende des sites de restriction si elle existe
  restriction_legend <- ""
  if (!is.null(restriction_sites_list) && length(restriction_sites_list) > 0) {
    restriction_legend <- generate_restriction_legend_formatted(restriction_sites_list)
  }

  # Retourner le contenu formaté pour le conteneur de légende
  return(paste0(legend_content, restriction_legend))
}

# ==================== GÉNÉRATION DES ALIGNEMENTS ====================

generate_colored_alignment <- function(pattern_seq, subject_seq, annotation_seq, start_position, colors, filename, restriction_positions = NULL) {
  line_length <- 60
  sequence_length <- nchar(subject_seq)

  html_lines <- character()
  text_lines <- character()

  html_lines <- c(html_lines, "<div class='alignment-block'>")
  html_lines <- c(html_lines, paste0("<div class='alignment-title'>Alignement avec ", filename, "</div>"))

  text_lines <- c(text_lines, paste0("=== Alignement avec ", filename, " ==="))

  for (start in seq(1, sequence_length, by = line_length)) {
    end <- min(start + line_length - 1, sequence_length)

    subject_segment <- substr(subject_seq, start, end)
    pattern_segment <- substr(pattern_seq, start, end)
    annotation_segment <- substr(annotation_seq, start, end)

    subject_start <- start
    subject_end <- end

    pattern_bases_before <- nchar(gsub("-", "", substr(pattern_seq, 1, start - 1)))
    pattern_bases_in_segment <- nchar(gsub("-", "", pattern_segment))
    pattern_start <- pattern_bases_before + 1
    pattern_end <- pattern_bases_before + pattern_bases_in_segment

    if (pattern_bases_in_segment == 0) {
      pattern_end <- pattern_bases_before
    }

    segment_colors <- if (!is.null(colors) && length(colors) >= end) {
      colors[start:end]
    } else {
      NULL
    }

    segment_restrictions <- if (!is.null(restriction_positions) && length(restriction_positions) >= end) {
      restriction_positions[start:end]
    } else {
      NULL
    }

    alignment_table <- create_colored_alignment_table(
      subject_segment, pattern_segment, annotation_segment,
      subject_start, subject_end, pattern_start, pattern_end,
      segment_colors, segment_restrictions
    )

    html_lines <- c(html_lines, alignment_table)

    text_lines <- c(text_lines, "")
    text_lines <- c(text_lines, sprintf("Seq_1 %6d %s %d", subject_start, subject_segment, subject_end))
    text_lines <- c(text_lines, sprintf("       %s", annotation_segment))
    text_lines <- c(text_lines, sprintf("Seq_2 %6d %s %d", pattern_start, pattern_segment, pattern_end))
  }

  html_lines <- c(html_lines, "</div>")

  return(list(
    html = paste(html_lines, collapse = "\n"),
    text = paste(text_lines, collapse = "\n")
  ))
}

create_colored_alignment_table <- function(subject_segment, pattern_segment, annotation_segment,
                                           subject_start, subject_end, pattern_start, pattern_end,
                                           segment_colors, restriction_positions = NULL) {
  subject_chars <- strsplit(subject_segment, "")[[1]]
  pattern_chars <- strsplit(pattern_segment, "")[[1]]
  annotation_chars <- strsplit(annotation_segment, "")[[1]]

  base_cell_style <- "font-family: Courier New, monospace; text-align: center; padding: 0; margin: 0; border: 0; width: 1ch; min-width: 1ch; max-width: 1ch;"

  subject_cells <- create_colored_sequence_cells(subject_chars, segment_colors, base_cell_style, restriction_positions)
  annotation_cells <- create_sequence_cells(annotation_chars, base_cell_style)
  pattern_cells <- create_sequence_cells(pattern_chars, base_cell_style)

  prefix_style <- "font-family: Courier New, monospace; padding: 0; margin: 0; text-align: left; width: 100px; min-width: 100px; max-width: 100px; white-space: nowrap;"
  suffix_style <- "font-family: Courier New, monospace; padding: 0; margin: 0; text-align: left; width: 60px; min-width: 60px; max-width: 60px;"

  table_style <- "border-collapse: collapse; margin: 10px 0; font-family: Courier New, monospace;"

  table_html <- paste0(
    '<table style="', table_style, '">',
    "<tr>",
    '<td style="', prefix_style, '">Seq_1 ', sprintf("%6d", subject_start), " </td>",
    paste(subject_cells, collapse = ""),
    '<td style="', suffix_style, '"> ', subject_end, "</td>",
    "</tr>",
    "<tr>",
    '<td style="', prefix_style, '"></td>',
    paste(annotation_cells, collapse = ""),
    '<td style="', suffix_style, '"></td>',
    "</tr>",
    "<tr>",
    '<td style="', prefix_style, '">Seq_2 ', sprintf("%6d", pattern_start), " </td>",
    paste(pattern_cells, collapse = ""),
    '<td style="', suffix_style, '"> ', pattern_end, "</td>",
    "</tr>",
    "</table>"
  )

  return(table_html)
}

# Créer les cellules colorées pour une séquence avec sites de restriction
create_colored_sequence_cells <- function(sequence_chars, colors, base_style, restriction_positions = NULL) {
  cells <- character()

  for (i in seq_along(sequence_chars)) {
    style <- base_style

    # Appliquer les couleurs normales (comme avant)
    if (!is.null(colors) && i <= length(colors) && !is.na(colors[i])) {
      style <- paste0(style, " color: ", colors[i], "; font-weight: bold;")
    }

    # Ajouter le fond gris si c'est un site de restriction
    if (!is.null(restriction_positions) && i <= length(restriction_positions) && restriction_positions[i]) {
      style <- paste0(style, " background-color: #D3D3D3; border: 1px solid #999; border-radius: 2px;")
    }

    cells <- c(cells, paste0("<td style='", style, "'>", sequence_chars[i], "</td>"))
  }
  return(cells)
}

# Créer les cellules normales pour une séquence
create_sequence_cells <- function(sequence_chars, base_style) {
  cells <- character()
  for (i in seq_along(sequence_chars)) {
    cells <- c(cells, paste0('<td style="', base_style, '">', sequence_chars[i], '</td>'))
  }
  return(cells)
}
