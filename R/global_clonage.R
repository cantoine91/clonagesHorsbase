# ==============================================================================
# GLOBAL_CLONAGE.R
# Fonctions utilitaires pour l'analyse de séquences et l'alignement
# Application Shiny HGX - Module Clonage - VERSION CORRIGÉE ESPACEMENT
# ==============================================================================

library(Biostrings)

# ==============================================================================
# FONCTIONS UTILITAIRES GÉNÉRALES
# ==============================================================================

#' Nettoyage et normalisation d'une séquence ADN
#' @param sequence_text Texte contenant la séquence à nettoyer
#' @return Séquence nettoyée en majuscules sans espaces ni retours à la ligne
clean_sequence <- function(sequence_text) {
  sequence_text <- gsub("[\r\n ]", "", sequence_text)
  toupper(sequence_text)
}

#' Validation de l'existence d'un répertoire
#' @param path Chemin du répertoire à valider
#' @return TRUE si le répertoire existe, FALSE sinon
validate_directory <- function(path) {
  return(dir.exists(path))
}

#' Nettoyage et validation des mots-clés de recherche
#' @param keyword Mot-clé à nettoyer
#' @return Mot-clé nettoyé sans espaces de début/fin
clean_search_keyword <- function(keyword) {
  if (is.null(keyword)) return("")
  trimws(keyword)
}

#' Formatage des résultats de recherche pour affichage
#' @param folders Vecteur des dossiers trouvés
#' @param files Vecteur des fichiers trouvés
#' @return Chaîne formatée décrivant les résultats
format_search_results <- function(folders, files) {
  if (length(folders) == 0) {
    return("Aucun résultat trouvé")
  }

  result <- paste("Dossiers trouvés:", length(folders), "\n")
  if (length(files) > 0) {
    result <- paste0(result, "Fichiers .seq trouvés:", length(files))
  } else {
    result <- paste0(result, "Aucun fichier .seq correspondant")
  }

  return(result)
}

#' Calcul de la région d'affichage basée sur les sites de restriction
#' @param restriction_sites_list Liste des sites de restriction
#' @param sequence_length Longueur totale de la séquence
#' @param context_length Nombre de nucléotides de contexte (défaut: 200)
#' @return Liste avec start et end de la région à afficher
calculate_restriction_display_region <- function(restriction_sites_list, sequence_length, context_length = 200) {
  if (length(restriction_sites_list) == 0) {
    return(list(start = 1, end = sequence_length))
  }

  # Collecter tous les sites de restriction
  all_sites <- c()
  for (enzyme_name in names(restriction_sites_list)) {
    sites <- restriction_sites_list[[enzyme_name]]
    if (length(sites) > 0) {
      all_sites <- c(all_sites, sites)
    }
  }

  if (length(all_sites) == 0) {
    return(list(start = 1, end = sequence_length))
  }

  # Trier les sites
  all_sites <- sort(all_sites)

  if (length(all_sites) == 1) {
    # Un seul site : centrer autour de ce site
    center <- all_sites[1]
    start_pos <- max(1, center - context_length)
    end_pos <- min(sequence_length, center + context_length)
  } else {
    # Plusieurs sites : du premier au dernier avec contexte
    first_site <- all_sites[1]
    last_site <- all_sites[length(all_sites)]

    start_pos <- max(1, first_site - context_length)
    end_pos <- min(sequence_length, last_site + context_length)
  }

  return(list(start = start_pos, end = end_pos))
}

# ==============================================================================
# ENZYMES DE RESTRICTION
# ==============================================================================

#' Base de données des enzymes de restriction communes
#' @return Liste nommée des séquences de reconnaissance des enzymes
get_restriction_enzymes <- function() {
  return(list(
    "AflII" = "CTTAAG",
    "ApaI" = "GGGCCC",
    "AscI" = "GGCGCGCC",
    "AvrII" = "CCTAGG",
    "BamHI" = "GGATCC",
    "BglII" = "AGATCT",
    "BspEI" = "TCCGGA",
    "BssHII" = "GCGCGC",
    "ClaI" = "ATCGAT",
    "EcoRI" = "GAATTC",
    "HindIII" = "AAGCTT",
    "HpaI" = "GTTAAC",
    "KpnI" = "GGTACC",
    "MluI" = "ACGCGT",
    "NcoI" = "CCATGG",
    "NdeI" = "CATATG",
    "NheI" = "GCTAGC",
    "NotI" = "GCGGCCGC",
    "PacI" = "TTAATTAA",
    "PmeI" = "GTTTAAAC",
    "PstI" = "CTGCAG",
    "SacI" = "GAGCTC",
    "SacII" = "CCGCGG",
    "SalI" = "GTCGAC",
    "SbfI" = "CCTGCAGG",
    "ScaI" = "AGTACT",
    "SfiI" = "GGCC.....GGCC", # Pattern spécial pour SfiI
    "SmaI" = "CCCGGG",
    "SpeI" = "ACTAGT",
    "SrfI" = "GCCCGGGC",
    "StuI" = "AGGCCT",
    "SwaI" = "ATTTAAAT",
    "XbaI" = "TCTAGA",
    "XhoI" = "CTCGAG",
    "XmaI" = "CCCGGG"
  ))
}

#' Recherche de sites de restriction dans une séquence
#' @param sequence Séquence ADN (DNAString ou caractère)
#' @param enzyme_sequence Séquence de reconnaissance de l'enzyme
#' @return Vecteur des positions de début des sites trouvés
find_restriction_sites <- function(sequence, enzyme_sequence) {
  if (is.null(sequence) || is.null(enzyme_sequence) || enzyme_sequence == "") {
    return(integer(0))
  }

  # Conversion en chaîne de caractères si nécessaire
  if (class(sequence)[1] == "DNAString") {
    sequence <- as.character(sequence)
  }

  sequence <- toupper(sequence)
  enzyme_sequence <- toupper(enzyme_sequence)

  # Gestion spéciale pour SfiI
  if (enzyme_sequence == "GGCC.....GGCC") {
    # Pattern regex pour SfiI : GGCC suivi de 5 n'importe quels nucléotides puis GGCC
    pattern <- "GGCC[ATCG]{5}GGCC"
    matches <- gregexpr(pattern, sequence, perl = TRUE)[[1]]

    if (matches[1] == -1) {
      return(integer(0))
    }
    return(as.integer(matches))
  }

  # Recherche normale pour les autres enzymes
  matches <- gregexpr(enzyme_sequence, sequence, fixed = TRUE)[[1]]

  if (matches[1] == -1) {
    return(integer(0))
  }

  return(as.integer(matches))
}


#' Construction d'une carte des positions de sites de restriction
#' @param sequence_length Longueur de la séquence
#' @param restriction_sites Liste des sites de restriction par enzyme
#' @param alignment_start Position de début de l'alignement (défaut: 1)
#' @return Vecteur logique indiquant les positions des sites de restriction
build_restriction_position_map <- function(sequence_length, restriction_sites, alignment_start = 1) {
  is_restriction <- rep(FALSE, sequence_length)

  if (!is.null(restriction_sites) && length(restriction_sites) > 0) {
    enzymes <- get_restriction_enzymes()

    for (enzyme_name in names(restriction_sites)) {
      sites <- restriction_sites[[enzyme_name]]
      if (length(sites) > 0) {

        # Gestion spéciale pour SfiI (site de 13 bp)
        if (enzyme_name == "SfiI") {
          site_length <- 13  # GGCC + 5 nucléotides + GGCC
        } else {
          site_length <- nchar(enzymes[[enzyme_name]])
        }

        for (site_pos in sites) {
          # Calcul des positions relatives à la fenêtre d'affichage
          for (i in 0:(site_length - 1)) {
            global_pos <- site_pos + i
            # Position relative dans la fenêtre d'affichage
            relative_pos <- global_pos - alignment_start + 1

            # Vérifier que la position est dans la fenêtre
            if (relative_pos > 0 && relative_pos <= sequence_length) {
              is_restriction[relative_pos] <- TRUE
            }
          }
        }
      }
    }
  }

  return(is_restriction)
}


#' Génération de la légende HTML pour les sites de restriction (mise à jour)
#' @param restriction_sites_list Liste des sites de restriction
#' @return HTML formaté pour la légende des sites de restriction
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

      # Affichage spécial pour SfiI
      if (enzyme_name == "SfiI") {
        enzyme_display <- "GGCCNNNNNGGCC"
      } else {
        enzyme_display <- enzyme_seq
      }

      sites_text <- paste(sites, collapse = ", ")
      legend_items <- c(legend_items,
                        paste0("<div style='margin-bottom: 5px;'><span style='color:", color,
                               "; font-weight: bold; font-size: 16px;'>■</span> ",
                               "<span style='font-size: 12px;'>", enzyme_name, " (", enzyme_display,
                               ") - Sites: ", sites_text, "</span></div>"))
      site_index <- site_index + 1
    }
  }

  if (length(legend_items) > 0) {
    return(paste0("<h4 style='margin-top: 15px; color: #b22222;'>Sites de restriction</h4>",
                  paste(legend_items, collapse = "")))
  }

  return("")
}

# ==============================================================================
# FONCTIONS D'ALIGNEMENT
# ==============================================================================

#' Création d'une séquence pattern complète avec gaps
#' @param aligned_pattern Séquence pattern alignée
#' @param start_position Position de début dans la séquence de référence
#' @param total_length Longueur totale de la séquence de référence
#' @return Séquence complète avec gaps aux positions non alignées
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

#' Création d'une annotation complète avec espaces
#' @param aligned_annotation Annotation alignée
#' @param start_position Position de début dans la séquence de référence
#' @param total_length Longueur totale de la séquence de référence
#' @return Annotation complète avec espaces aux positions non alignées
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

#' Annotation des mutations entre deux séquences alignées
#' @param pattern_sequence Séquence pattern alignée
#' @param subject_sequence Séquence sujet alignée
#' @return Chaîne d'annotation avec "|" pour les correspondances exactes
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
    if (pattern_chars[i] == subject_chars[i] &&
        pattern_chars[i] != "-" && subject_chars[i] != "-") {
      annotation[i] <- "|"
    } else {
      annotation[i] <- " "
    }
  }

  return(paste(annotation, collapse = ""))
}


#' Création des cellules d'annotation avec surlignage des mutations (VERSION CORRIGÉE)
#' @param annotation_chars Caractères d'annotation
#' @param base_style Style CSS de base
#' @param start_position Position de début
#' @param sequence_type Type de séquence
#' @return Vecteur de cellules HTML formatées avec surlignage des mutations
create_annotation_cells_with_mutation_highlight <- function(annotation_chars, base_style, start_position = 1, sequence_type = "annotation") {
  cells <- character()

  # Trouver la première et dernière position avec '|'
  match_positions <- which(annotation_chars == "|")

  if (length(match_positions) > 0) {
    first_match <- match_positions[1]
    last_match <- match_positions[length(match_positions)]
  } else {
    # Si pas de '|', pas de surlignage
    first_match <- 0
    last_match <- 0
  }

  for (i in seq_along(annotation_chars)) {
    style <- base_style

    # Ajouter le fond rouge si :
    # 1. On est dans la zone d'alignement (du premier au dernier '|')
    # 2. ET le caractère n'est pas '|' (donc c'est une mutation)
    if (first_match > 0 &&
        i >= first_match &&
        i <= last_match &&
        annotation_chars[i] != "|") {
      # Fond rouge pour les mutations
      style <- paste0(style, " background-color: #ffcccc; border: 1px solid #ff9999; border-radius: 2px;")
    }

    # Créer la cellule
    cells <- c(cells, paste0('<td style="', style, '">', annotation_chars[i], '</td>'))
  }

  return(cells)
}

#' Création d'un tableau HTML coloré pour un segment d'alignement (VERSION CORRIGÉE)
#' @param subject_segment Segment de la séquence sujet
#' @param pattern_segment Segment de la séquence pattern
#' @param annotation_segment Segment d'annotation
#' @param subject_start Position de début du sujet
#' @param subject_end Position de fin du sujet
#' @param pattern_start Position de début du pattern
#' @param pattern_end Position de fin du pattern
#' @param segment_colors Couleurs pour ce segment
#' @param restriction_positions Positions des restrictions pour ce segment
#' @return HTML du tableau formaté
create_colored_alignment_table <- function(subject_segment, pattern_segment, annotation_segment,
                                           subject_start, subject_end, pattern_start, pattern_end,
                                           segment_colors, restriction_positions = NULL) {
  # Décomposition en caractères
  subject_chars <- strsplit(subject_segment, "")[[1]]
  pattern_chars <- strsplit(pattern_segment, "")[[1]]
  annotation_chars <- strsplit(annotation_segment, "")[[1]]

  # Styles de base
  base_cell_style <- "font-family: Courier New, monospace; text-align: center; padding: 0; margin: 0; border: 0; width: 1ch; min-width: 1ch; max-width: 1ch;"
  prefix_style <- "font-family: Courier New, monospace; padding: 0; margin: 0; text-align: left; width: 100px; min-width: 100px; max-width: 100px; white-space: nowrap;"
  suffix_style <- "font-family: Courier New, monospace; padding: 0 0 0 8px; margin: 0; text-align: left; width: 60px; min-width: 60px; max-width: 60px;"
  table_style <- "border-collapse: collapse; margin: 10px 0; font-family: Courier New, monospace;"

  # Création des cellules
  subject_cells <- create_colored_sequence_cells(subject_chars, segment_colors, base_cell_style,
                                                 restriction_positions, subject_start, "subject")

  # UTILISATION DE LA FONCTION CORRIGÉE pour les annotations
  annotation_cells <- create_annotation_cells_with_mutation_highlight(annotation_chars, base_cell_style, 1, "annotation")

  pattern_cells <- create_sequence_cells(pattern_chars, base_cell_style, pattern_start, "pattern")

  # Assemblage du tableau HTML
  table_html <- paste0(
    '<table style="', table_style, '">',
    "<tr>",
    '<td style="', prefix_style, '">Seq_1 ', sprintf("%6d", subject_start), " </td>",
    paste(subject_cells, collapse = ""),
    '<td style="', suffix_style, '">&nbsp;', subject_end, "</td>",
    "</tr>",
    "<tr>",
    '<td style="', prefix_style, '"></td>',
    paste(annotation_cells, collapse = ""),
    '<td style="', suffix_style, '"></td>',
    "</tr>",
    "<tr>",
    '<td style="', prefix_style, '">Seq_2 ', sprintf("%6d", pattern_start), " </td>",
    paste(pattern_cells, collapse = ""),
    '<td style="', suffix_style, '">&nbsp;', pattern_end, "</td>",
    "</tr>",
    "</table>"
  )

  return(table_html)
}

# ==============================================================================
# PARSING DES FEATURES GENBANK
# ==============================================================================

#' Extraction de la couleur hexadécimale des features SerialCloner
#' @param feature_lines Lignes de la feature GenBank
#' @return Couleur hexadécimale (défaut: "#000000")
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

#' Parsing des features GenBank de manière robuste
#' @param features_lines Lignes contenant les features du fichier GenBank
#' @return Liste des features parsées avec type, position, couleur et nom
parse_genbank_features <- function(features_lines) {
  features_list <- list()
  current_feature <- NULL
  current_feature_lines <- character()

  for (i in seq_along(features_lines)) {
    line <- features_lines[i]
    is_new_feature_with_type <- grepl("^\\s{5,}[a-zA-Z_]+\\s+", line)
    is_position_only <- grepl("^\\s+\\d+\\.\\.\\d+", line)

    if (is_new_feature_with_type || is_position_only) {
      # Sauvegarder la feature précédente
      if (!is.null(current_feature)) {
        current_feature$color <- extract_serialcloner_color(current_feature_lines)
        features_list[[length(features_list) + 1]] <- current_feature
      }

      # Créer une nouvelle feature
      feature_type <- if (is_new_feature_with_type) extract_feature_type(line) else "misc_feature"
      position_raw <- extract_feature_position(line)

      current_feature <- list(
        type = feature_type,
        position_raw = position_raw,
        color = "#000000",
        name = ""
      )
      current_feature_lines <- character()

      # Gestion des positions orphelines
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

    # Accumulation des lignes et extraction du nom
    if (!is.null(current_feature)) {
      current_feature_lines <- c(current_feature_lines, line)
      feature_name <- extract_feature_name(line)
      if (!is.null(feature_name) && current_feature$name == "") {
        current_feature$name <- feature_name
      }
    }
  }

  # Sauvegarder la dernière feature
  if (!is.null(current_feature)) {
    current_feature$color <- extract_serialcloner_color(current_feature_lines)
    features_list[[length(features_list) + 1]] <- current_feature
  }

  return(features_list)
}

#' Extraction du type de feature
#' @param line Ligne contenant la déclaration de feature
#' @return Type de feature ou "unknown"
extract_feature_type <- function(line) {
  type_match <- regmatches(line, regexpr("^\\s+([a-zA-Z_]+)", line))
  if (length(type_match) > 0) {
    return(gsub("^\\s+", "", type_match))
  }
  return("unknown")
}

#' Extraction de la position de la feature
#' @param line Ligne contenant la position
#' @return Position au format "start..end" ou NA
extract_feature_position <- function(line) {
  # Position normale
  position_match <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
  if (length(position_match) > 0) {
    return(position_match)
  }

  # Position complement
  complement_match <- regmatches(line, regexpr("complement\\(\\d+\\.\\.\\d+\\)", line))
  if (length(complement_match) > 0) {
    position_match <- regmatches(complement_match, regexpr("\\d+\\.\\.\\d+", complement_match))
    if (length(position_match) > 0) {
      return(position_match)
    }
  }

  return(NA)
}

#' Extraction du nom de la feature (label ou gene)
#' @param line Ligne contenant le label ou gene
#' @return Nom de la feature ou NULL
extract_feature_name <- function(line) {
  # Gestion des labels
  if (grepl("/label\\s*=", line)) {
    if (grepl('/label\\s*=\\s*"[^"]*"', line)) {
      label_match <- regmatches(line, regexpr('/label\\s*=\\s*"[^"]*"', line))
      return(gsub('/label\\s*=\\s*"([^"]*)"', '\\1', label_match))
    } else {
      label_match <- regmatches(line, regexpr('/label\\s*=\\s*.+', line))
      if (length(label_match) > 0) {
        return(gsub('/label\\s*=\\s*', '', label_match))
      }
    }
  }

  # Gestion des genes
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

# ==============================================================================
# GESTION DES COULEURS
# ==============================================================================

#' Configuration des couleurs par défaut pour les features connues
#' @return Liste nommée des couleurs par défaut
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

#' Obtention de la couleur par nom de feature
#' @param feature_name Nom de la feature
#' @return Couleur hexadécimale (défaut: "#888888")
get_color_by_feature_name <- function(feature_name) {
  default_colors <- get_default_feature_colors()
  color <- default_colors[[feature_name]]
  if (is.null(color)) {
    return("#888888")
  }
  return(color)
}

#' Construction du vecteur de couleurs pour la séquence
#' @param features_lines Lignes des features GenBank
#' @param alignment_start Position de début de l'alignement
#' @param sequence_length Longueur de la séquence
#' @return Vecteur de couleurs pour chaque position
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

    # Détermination de la couleur
    color_to_use <- feature$color
    if (color_to_use == "#000000") {
      color_to_use <- get_color_by_feature_name(feature$name)
    }

    # Application de la couleur sur la plage
    for (i in seq_len(sequence_length)) {
      reference_position <- alignment_start + i - 1
      if (reference_position >= feature_start && reference_position <= feature_end) {
        color_map[i] <- color_to_use
      }
    }
  }

  return(color_map)
}

#' Génération de la légende HTML des couleurs
#' @param features_lines Lignes des features GenBank
#' @param restriction_sites_list Liste des sites de restriction (optionnel)
#' @return HTML formaté pour la légende complète
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
                      paste0("<div style='margin-bottom: 5px;'><span style='color:", color_to_use,
                             "; font-weight: bold; font-size: 16px;'>■</span> ",
                             "<span style='font-size: 12px;'>", display_name, " (",
                             position_bounds[1], "-", position_bounds[2], ")</span></div>"))
  }

  # Construction de la légende des features
  legend_content <- if (length(legend_items) > 0) {
    paste0("<h4 style='margin-top: 0; color: #b22222;'>Légende des couleurs</h4>",
           paste(legend_items, collapse = ""))
  } else {
    "<h4 style='margin-top: 0; color: #b22222;'>Aucune annotation trouvée</h4>"
  }

  # Ajout de la légende des sites de restriction
  restriction_legend <- ""
  if (!is.null(restriction_sites_list) && length(restriction_sites_list) > 0) {
    restriction_legend <- generate_restriction_legend_formatted(restriction_sites_list)
  }

  return(paste0(legend_content, restriction_legend))
}

# ==============================================================================
# GÉNÉRATION DES ALIGNEMENTS HTML
# ==============================================================================

#' Génération d'un alignement coloré au format HTML et texte
#' @param pattern_seq Séquence pattern complète
#' @param subject_seq Séquence sujet complète
#' @param annotation_seq Annotation d'alignement complète
#' @param start_position Position de début de l'alignement
#' @param colors Vecteur de couleurs pour chaque position
#' @param filename Nom du fichier pour l'en-tête
#' @param restriction_positions Positions des sites de restriction (optionnel)
#' @return Liste avec versions HTML et texte de l'alignement
generate_colored_alignment <- function(pattern_seq, subject_seq, annotation_seq, start_position,
                                       colors, filename, restriction_positions = NULL) {
  line_length <- 60
  sequence_length <- nchar(subject_seq)

  html_lines <- character()
  text_lines <- character()

  # En-têtes
  html_lines <- c(html_lines, "<div class='alignment-block'>")
  html_lines <- c(html_lines, paste0("<div class='alignment-title'>Alignement avec ", filename, "</div>"))
  text_lines <- c(text_lines, paste0("=== Alignement avec ", filename, " ==="))

  # Génération par blocs de 60 caractères
  for (start in seq(1, sequence_length, by = line_length)) {
    end <- min(start + line_length - 1, sequence_length)

    # Extraction des segments
    subject_segment <- substr(subject_seq, start, end)
    pattern_segment <- substr(pattern_seq, start, end)
    annotation_segment <- substr(annotation_seq, start, end)

    # Calcul des positions d'affichage pour le sujet (séquence de référence)
    subject_start <- start
    subject_end <- end

    # CORRECTION: Calcul correct des positions pour le pattern
    # Compter les nucléotides réels (non-gap) avant ce segment
    pattern_before_segment <- substr(pattern_seq, 1, start - 1)
    pattern_bases_before <- nchar(gsub("-", "", pattern_before_segment))

    # Compter les nucléotides réels dans ce segment
    pattern_bases_in_segment <- nchar(gsub("-", "", pattern_segment))

    # Calculer les positions de début et fin
    if (pattern_bases_in_segment == 0) {
      # Segment entièrement composé de gaps
      pattern_start <- pattern_bases_before
      pattern_end <- pattern_bases_before
    } else {
      pattern_start <- pattern_bases_before + 1
      pattern_end <- pattern_bases_before + pattern_bases_in_segment
    }

    # Extraction des couleurs et restrictions pour ce segment
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

    # Génération du tableau HTML pour ce segment
    alignment_table <- create_colored_alignment_table(
      subject_segment, pattern_segment, annotation_segment,
      subject_start, subject_end, pattern_start, pattern_end,
      segment_colors, segment_restrictions
    )

    html_lines <- c(html_lines, alignment_table)

    # Version texte - CORRECTION ESPACEMENT
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

#' Création d'un tableau HTML coloré pour un segment d'alignement
#' @param subject_segment Segment de la séquence sujet
#' @param pattern_segment Segment de la séquence pattern
#' @param annotation_segment Segment d'annotation
#' @param subject_start Position de début du sujet
#' @param subject_end Position de fin du sujet
#' @param pattern_start Position de début du pattern
#' @param pattern_end Position de fin du pattern
#' @param segment_colors Couleurs pour ce segment
#' @param restriction_positions Positions des restrictions pour ce segment
#' @return HTML du tableau formaté
create_colored_alignment_table <- function(subject_segment, pattern_segment, annotation_segment,
                                           subject_start, subject_end, pattern_start, pattern_end,
                                           segment_colors, restriction_positions = NULL) {
  # Décomposition en caractères
  subject_chars <- strsplit(subject_segment, "")[[1]]
  pattern_chars <- strsplit(pattern_segment, "")[[1]]
  annotation_chars <- strsplit(annotation_segment, "")[[1]]

  # Styles de base
  base_cell_style <- "font-family: Courier New, monospace; text-align: center; padding: 0; margin: 0; border: 0; width: 1ch; min-width: 1ch; max-width: 1ch;"
  prefix_style <- "font-family: Courier New, monospace; padding: 0; margin: 0; text-align: left; width: 100px; min-width: 100px; max-width: 100px; white-space: nowrap;"
  suffix_style <- "font-family: Courier New, monospace; padding: 0 0 0 8px; margin: 0; text-align: left; width: 60px; min-width: 60px; max-width: 60px;"
  table_style <- "border-collapse: collapse; margin: 10px 0; font-family: Courier New, monospace;"

  # Création des cellules
  subject_cells <- create_colored_sequence_cells(subject_chars, segment_colors, base_cell_style,
                                                 restriction_positions, subject_start, "subject")
  annotation_cells <- create_sequence_cells(annotation_chars, base_cell_style, 1, "annotation")
  pattern_cells <- create_sequence_cells(pattern_chars, base_cell_style, pattern_start, "pattern")

  # Assemblage du tableau HTML - CORRECTION ESPACEMENT FORCÉ
  table_html <- paste0(
    '<table style="', table_style, '">',
    "<tr>",
    '<td style="', prefix_style, '">Seq_1 ', sprintf("%6d", subject_start), " </td>",
    paste(subject_cells, collapse = ""),
    '<td style="', suffix_style, '">&nbsp;', subject_end, "</td>",
    "</tr>",
    "<tr>",
    '<td style="', prefix_style, '"></td>',
    paste(annotation_cells, collapse = ""),
    '<td style="', suffix_style, '"></td>',
    "</tr>",
    "<tr>",
    '<td style="', prefix_style, '">Seq_2 ', sprintf("%6d", pattern_start), " </td>",
    paste(pattern_cells, collapse = ""),
    '<td style="', suffix_style, '">&nbsp;', pattern_end, "</td>",
    "</tr>",
    "</table>"
  )

  return(table_html)
}

#' Création des cellules colorées avec tooltips pour la séquence sujet
#' @param sequence_chars Caractères de la séquence
#' @param colors Couleurs à appliquer
#' @param base_style Style CSS de base
#' @param restriction_positions Positions des sites de restriction
#' @param start_position Position de début dans la séquence
#' @param sequence_type Type de séquence ("subject" ou "pattern")
#' @return Vecteur de cellules HTML formatées
create_colored_sequence_cells <- function(sequence_chars, colors, base_style, restriction_positions = NULL,
                                          start_position = 1, sequence_type = "subject") {
  cells <- character()

  for (i in seq_along(sequence_chars)) {
    style <- paste0(base_style, " position: relative; cursor: pointer;")

    # Calcul de la position réelle dans la séquence
    real_position <- start_position + i - 1

    # Application des couleurs normales
    if (!is.null(colors) && i <= length(colors) && !is.na(colors[i])) {
      style <- paste0(style, " color: ", colors[i], "; font-weight: bold;")
    }

    # Ajout du fond gris pour les sites de restriction
    if (!is.null(restriction_positions) && i <= length(restriction_positions) && restriction_positions[i]) {
      style <- paste0(style, " background-color: #D3D3D3; border: 1px solid #999; border-radius: 2px;")
    }

    # Création du tooltip avec informations de position
    tooltip_text <- paste0("Position: ", real_position, " | Nucléotide: ", sequence_chars[i])
    if (sequence_type == "pattern") {
      tooltip_text <- paste0(tooltip_text, " (Séquence test)")
    } else {
      tooltip_text <- paste0(tooltip_text, " (Séquence référence)")
    }

    tooltip_html <- paste0('<span class="tooltip">', tooltip_text, '</span>')

    cells <- c(cells, paste0('<td class="nucleotide-cell" style="', style, '">',
                             sequence_chars[i], tooltip_html, '</td>'))
  }
  return(cells)
}

#' Version alternative simplifiée pour créer les cellules d'annotation
create_sequence_cells <- function(sequence_chars, base_style, start_position = 1, sequence_type = "pattern") {
  cells <- character()

  # Si c'est une annotation, appliquer le surlignage des mutations
  if (sequence_type == "annotation") {
    # Trouver la première et dernière position avec '|'
    match_positions <- which(sequence_chars == "|")

    if (length(match_positions) > 0) {
      first_match <- match_positions[1]
      last_match <- match_positions[length(match_positions)]
    } else {
      first_match <- 0
      last_match <- 0
    }

    for (i in seq_along(sequence_chars)) {
      style <- paste0(base_style, " position: relative; cursor: pointer;")

      # Ajouter le fond rouge pour les mutations dans la zone d'alignement
      if (first_match > 0 &&
          i >= first_match &&
          i <= last_match &&
          sequence_chars[i] != "|") {
        style <- paste0(style, " background-color: #ffcccc; border: 1px solid #ff9999; border-radius: 2px;")
      }

      cells <- c(cells, paste0('<td style="', style, '">', sequence_chars[i], '</td>'))
    }

    return(cells)
  }

  # Pour les séquences pattern normales
  for (i in seq_along(sequence_chars)) {
    style <- paste0(base_style, " position: relative; cursor: pointer;")

    # Gestion spécifique pour les séquences pattern
    if (sequence_type == "pattern") {
      # Calcul de position en excluant les gaps
      if (i == 1) {
        if (sequence_chars[i] == "-") {
          real_position <- NA
        } else {
          real_position <- start_position
        }
      } else {
        chars_before <- sequence_chars[1:(i-1)]
        non_gap_before <- sum(chars_before != "-")
        if (sequence_chars[i] == "-") {
          real_position <- NA
        } else {
          real_position <- start_position + non_gap_before
        }
      }

      if (is.na(real_position) || sequence_chars[i] == "-") {
        tooltip_text <- "Gap (insertion dans la référence)"
      } else {
        tooltip_text <- paste0("Position: ", real_position, " | Nucléotide: ", sequence_chars[i], " (Séquence test)")
      }
    } else {
      tooltip_text <- ""
    }

    # Création du tooltip si nécessaire
    tooltip_html <- if (tooltip_text != "" && sequence_chars[i] != " ") {
      paste0('<span class="tooltip">', tooltip_text, '</span>')
    } else {
      ""
    }

    cell_class <- if (tooltip_text != "" && sequence_chars[i] != " ") "nucleotide-cell" else ""

    cells <- c(cells, paste0('<td class="', cell_class, '" style="', style, '">',
                             sequence_chars[i], tooltip_html, '</td>'))
  }
  return(cells)
}

# ==============================================================================
# Gestion des groupes de clones
# ==============================================================================

#' Extraction du groupe de clone par avant-dernier et dernier underscore
#' @param filename Nom du fichier .seq
#' @return Liste avec prefix, group et full_group
extract_clone_group <- function(filename) {
  # Enlever l'extension .seq
  base_name <- gsub("\\.seq$", "", basename(filename))

  # Séparer par underscore
  parts <- strsplit(base_name, "_")[[1]]

  if (length(parts) >= 2) {
    # Ce qu'il y a entre l'avant-dernier et le dernier underscore
    clone_id <- parts[length(parts) - 1]  # Avant-dernier élément

    # Utiliser le clone_id comme groupe
    group_key <- clone_id

    # Le préfixe = tout sauf les deux derniers éléments
    if (length(parts) > 2) {
      prefix <- paste(parts[1:(length(parts)-2)], collapse = "_")
    } else {
      prefix <- parts[1]
    }

    return(list(
      prefix = prefix,
      group = group_key,
      full_group = paste(prefix, group_key, sep = "_")
    ))
  } else {
    # Si pas assez d'underscores, utiliser le nom complet
    return(list(
      prefix = base_name,
      group = base_name,
      full_group = base_name
    ))
  }
}

#' Organisation des fichiers par groupes - VERSION SIMPLIFIÉE
#' @param file_paths Vecteur des chemins complets des fichiers
#' @param file_names Vecteur des noms d'affichage
#' @return Liste organisée par groupes
organize_files_by_groups <- function(file_paths, file_names) {
  if (length(file_paths) == 0) {
    return(list())
  }

  groups <- list()

  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    display_name <- file_names[i]

    # Extraire le groupe
    clone_info <- extract_clone_group(file_path)
    group_key <- clone_info$group  # Ce sera "A-1", "A-2", etc.

    if (is.null(groups[[group_key]])) {
      groups[[group_key]] <- list(
        files = character(),
        paths = character(),
        display_names = character()
      )
    }

    groups[[group_key]]$files <- c(groups[[group_key]]$files, basename(file_path))
    groups[[group_key]]$paths <- c(groups[[group_key]]$paths, file_path)
    groups[[group_key]]$display_names <- c(groups[[group_key]]$display_names, display_name)
  }

  return(groups)
}
