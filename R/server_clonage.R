# ==============================================================================
# SERVER_CLONAGE.R
# Logique serveur pour l'application Shiny HGX - Module Clonage
# VERSION COMPLÈTE AVEC TOUTES LES FONCTIONS - ALIGNEMENT PDF CORRIGÉ
# ==============================================================================

library(shiny)
library(Biostrings)

server_clonage <- function(input, output, session) {

  # ==============================================================================
  # VARIABLES RÉACTIVES
  # ==============================================================================

  data_xdna <- reactiveValues(
    seq = NULL,
    features = NULL
  )

  alignment_data <- reactiveValues(
    results = NULL,
    text_version = NULL
  )

  search_data <- reactiveValues(
    folders_found = character(),
    seq_files_found = character(),
    seq_files_paths = character(),
    search_in_progress = FALSE,
    stop_search = FALSE
  )

  alignment_progress <- reactiveValues(
    in_progress = FALSE,
    complete = FALSE,
    status_message = "Prêt"
  )

  # ==============================================================================
  # CONFIGURATION DES CHEMINS
  # ==============================================================================

  if (dir.exists("/data/production/SEQ")) {
    xdna_dir <- "/data/production/SEQ/Atest_cae"
    seq_base_dir <- "/data/production/SEQ"
  } else if (dir.exists("../data/production/SEQ")) {
    xdna_dir <- "../data/production/SEQ/Atest_cae"
    seq_base_dir <- "../data/production/SEQ"
  } else {
    xdna_dir <- "P:/SEQ/Atest_cae"
    seq_base_dir <- "P:/SEQ"
  }

  cat("📁 Chemins configurés:\n")
  cat("   - GenBank (.gb):", xdna_dir, "\n")
  cat("   - Séquences (.seq):", seq_base_dir, "\n")

  # ==============================================================================
  # FONCTIONS DE RECHERCHE
  # ==============================================================================

  search_ultra_fast <- function(plate_keyword, base_dir = "P:/SEQ") {
    if (is.null(plate_keyword) || plate_keyword == "") {
      return(character())
    }

    if (!dir.exists(base_dir)) {
      return(character())
    }

    withProgress(message = 'Recherche ultra-rapide...', value = 0, {
      tryCatch({
        plate_keyword <- trimws(plate_keyword)
        incProgress(0.3, detail = paste("Recherche directe de '*", plate_keyword, "*'...", sep=""))

        search_pattern <- paste0("*", plate_keyword, "*")
        matching_folders <- Sys.glob(file.path(base_dir, search_pattern))
        matching_folders <- matching_folders[file.info(matching_folders)$isdir]

        if (length(matching_folders) == 0) {
          incProgress(1, detail = "Aucun dossier trouvé")
          return(character())
        }

        incProgress(0.3, detail = paste(length(matching_folders), "dossiers trouvés"))
        matching_folders <- sort(matching_folders, decreasing = TRUE)
        selected_folder <- matching_folders[1]

        incProgress(0.4, detail = paste("Sélectionné:", basename(selected_folder)))
        return(selected_folder)

      }, error = function(e) {
        incProgress(1, detail = paste("Erreur:", e$message))
        return(character())
      })
    })
  }

  search_all_seq_in_folder <- function(folder_path, seq_keyword) {
    if (length(folder_path) == 0 || is.null(seq_keyword) || seq_keyword == "") {
      return(list(files = character(), paths = character()))
    }

    tryCatch({
      seq_keyword <- trimws(seq_keyword)

      if (!dir.exists(folder_path)) {
        return(list(files = character(), paths = character()))
      }

      seq_files <- list.files(folder_path, pattern = "\\.seq$",
                              recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

      if (length(seq_files) == 0) {
        return(list(files = character(), paths = character()))
      }

      file_names <- basename(seq_files)
      matching_indices <- grep(seq_keyword, file_names, ignore.case = TRUE)

      if (length(matching_indices) == 0) {
        return(list(files = character(), paths = character()))
      }

      matching_files <- seq_files[matching_indices]
      folder_name <- basename(folder_path)
      display_names <- character()

      for (i in seq_along(matching_files)) {
        relative_path <- gsub(paste0("^", gsub("([\\(\\)\\[\\]\\{\\}\\^\\$\\*\\+\\?\\|\\\\])",
                                               "\\\\\\1", folder_path), "[\\/\\\\]?"), "", matching_files[i])

        if (dirname(relative_path) != ".") {
          display_names[i] <- paste0("[", folder_name, "] ", dirname(relative_path), " → ", basename(relative_path))
        } else {
          display_names[i] <- paste0("[", folder_name, "] ", basename(relative_path))
        }
      }

      return(list(files = display_names, paths = matching_files))

    }, error = function(e) {
      return(list(files = character(), paths = character()))
    })
  }

  # ==============================================================================
  # GESTION DES FICHIERS GENBANK
  # ==============================================================================

  get_available_gb_files <- function() {
    gb_files <- list.files(xdna_dir, pattern = "\\.gb$", full.names = FALSE)
    return(gb_files)
  }

  observe({
    gb_files <- get_available_gb_files()
    updateSelectInput(session, "carte_xdna", choices = gb_files)
  })

  observeEvent(input$refresh_files, {
    gb_files <- get_available_gb_files()
    updateSelectInput(session, "carte_xdna", choices = gb_files)
    showNotification("📁 Liste des fichiers GenBank mise à jour !", type = "message", duration = 2)
  })

  # ==============================================================================
  # SITES DE RESTRICTION
  # ==============================================================================

  restriction_sites <- reactive({
    req(data_xdna$seq)
    sites_list <- list()
    enzymes <- get_restriction_enzymes()

    if (!is.null(input$enzyme1) && input$enzyme1 != "") {
      sites1 <- find_restriction_sites(data_xdna$seq, enzymes[[input$enzyme1]])
      if (length(sites1) > 0) sites_list[[input$enzyme1]] <- sites1
    }

    if (!is.null(input$enzyme2) && input$enzyme2 != "") {
      sites2 <- find_restriction_sites(data_xdna$seq, enzymes[[input$enzyme2]])
      if (length(sites2) > 0) sites_list[[input$enzyme2]] <- sites2
    }

    return(sites_list)
  })

  output$restriction_info <- renderText({
    sites <- restriction_sites()
    if (length(sites) == 0) return("Aucun site trouvé")

    info <- sapply(names(sites), function(e) paste0(e, ": ", length(sites[[e]]), " site(s)"))
    paste(info, collapse = " | ")
  })

  # ==============================================================================
  # GESTION DE LA RECHERCHE
  # ==============================================================================

  observeEvent(input$stop_search, {
    search_data$stop_search <- TRUE
    search_data$search_in_progress <- FALSE
    showNotification("🛑 Recherche interrompue", type = "warning", duration = 2)
  })

  observeEvent(input$search_seq_btn, {
    req(input$plate_keyword, input$seq_keyword)

    if (nchar(trimws(input$plate_keyword)) == 0) {
      showNotification("⚠️ Veuillez saisir un nom de plaque", type = "warning", duration = 3)
      return()
    }

    if (nchar(trimws(input$seq_keyword)) == 0) {
      showNotification("⚠️ Veuillez saisir un mot-clé pour les fichiers .seq", type = "warning", duration = 3)
      return()
    }

    search_data$stop_search <- FALSE
    search_data$search_in_progress <- TRUE
    search_data$folders_found <- character()
    search_data$seq_files_found <- character()
    search_data$seq_files_paths <- character()
    updateSelectInput(session, "seq_files", choices = NULL)

    tryCatch({
      plate_folder <- search_ultra_fast(input$plate_keyword, seq_base_dir)

      if (length(plate_folder) == 0) {
        search_data$search_in_progress <- FALSE
        showNotification("❌ Aucun dossier trouvé pour cette plaque", type = "error", duration = 5)
        return()
      }

      search_data$folders_found <- plate_folder
      showNotification(paste("📂 Dossier sélectionné:", basename(plate_folder)), type = "message", duration = 3)

      if (!search_data$stop_search) {
        withProgress(message = 'Recherche fichiers .seq...', value = 0, {
          incProgress(0.5, detail = "Scan du dossier...")

          seq_results <- search_all_seq_in_folder(plate_folder, input$seq_keyword)
          search_data$seq_files_found <- seq_results$files
          search_data$seq_files_paths <- seq_results$paths

          incProgress(0.5, detail = paste(length(seq_results$files), "fichiers trouvés"))
          search_data$search_in_progress <- FALSE

          if (length(seq_results$files) == 0) {
            showNotification("❌ Aucun fichier .seq trouvé avec ce mot-clé", type = "warning", duration = 5)
          } else {
            choices_list <- setNames(seq_results$paths, seq_results$files)
            updateSelectInput(session, "seq_files", choices = choices_list)
            showNotification(paste("✅", length(seq_results$files), "fichier(s) .seq trouvé(s)"),
                             type = "message", duration = 3)
          }
        })
      }

    }, error = function(e) {
      search_data$search_in_progress <- FALSE
      showNotification(paste("❌ Erreur lors de la recherche:", as.character(e$message)),
                       type = "error", duration = 5)
    })
  })

  # ==============================================================================
  # OUTPUTS INTERFACE
  # ==============================================================================

  output$search_in_progress <- reactive({
    search_data$search_in_progress
  })
  outputOptions(output, "search_in_progress", suspendWhenHidden = FALSE)

  output$search_results <- renderText({
    req(length(search_data$folders_found) > 0 || length(search_data$seq_files_found) > 0)

    if (length(search_data$folders_found) > 0) {
      folder_name <- basename(search_data$folders_found[1])
      folders_text <- paste0("📂 Dossier sélectionné: ", folder_name)

      if (length(search_data$seq_files_found) > 0) {
        files_preview <- if (length(search_data$seq_files_found) > 5) {
          c(search_data$seq_files_found[1:5], "...")
        } else {
          search_data$seq_files_found
        }

        files_text <- paste0("<br>📄 Fichiers .seq trouvés (", length(search_data$seq_files_found), "): ",
                             paste(files_preview, collapse = ", "))
        return(paste0(folders_text, files_text))
      } else {
        return(paste0(folders_text, "<br>❌ Aucun fichier .seq trouvé avec le mot-clé '", input$seq_keyword, "'"))
      }
    } else {
      return("")
    }
  })

  output$seq_files_found <- reactive({
    length(search_data$seq_files_found) > 0
  })
  outputOptions(output, "seq_files_found", suspendWhenHidden = FALSE)

  # ==============================================================================
  # CHARGEMENT GENBANK
  # ==============================================================================

  observeEvent(input$align_btn, {
    req(input$carte_xdna)

    alignment_progress$in_progress <- TRUE
    alignment_progress$complete <- FALSE
    alignment_progress$status_message <- "Chargement du fichier GenBank..."

    Sys.sleep(0.1)

    fichier <- file.path(xdna_dir, input$carte_xdna)

    tryCatch({
      alignment_progress$status_message <- "Lecture du fichier GenBank..."
      gb_lines <- readLines(fichier, warn = FALSE, encoding = "UTF-8")
    }, error = function(e) {
      tryCatch({
        gb_lines <- readLines(fichier, warn = FALSE, encoding = "latin1")
      }, error = function(e2) {
        gb_lines <- readLines(fichier, warn = FALSE)
      })
    })

    alignment_progress$status_message <- "Nettoyage des données..."
    gb_lines <- iconv(gb_lines, to = "UTF-8", sub = "")
    gb_lines <- gb_lines[!is.na(gb_lines)]
    gb_lines <- gsub("[^\x01-\x7F]", "", gb_lines)

    alignment_progress$status_message <- "Extraction de la séquence..."
    origin_line <- tryCatch({
      grep("^ORIGIN", gb_lines, ignore.case = TRUE)
    }, warning = function(w) {
      which(grepl("^ORIGIN", gb_lines, ignore.case = TRUE))
    })

    if (length(origin_line) == 0) {
      alignment_progress$in_progress <- FALSE
      alignment_progress$status_message <- "Erreur: Section ORIGIN non trouvée"
      showNotification("Erreur: Section ORIGIN non trouvée dans le fichier GenBank",
                       type = "error", duration = 5)
      return()
    }

    alignment_progress$status_message <- "Extraction des annotations..."
    features_block <- gb_lines[1:(origin_line[1] - 1)]
    features_lines <- tryCatch({
      features_block[grep("^\\s{5}|^\\s{21}", features_block)]
    }, warning = function(w) {
      features_block[grepl("^\\s{5}|^\\s{21}", features_block)]
    })

    data_xdna$features <- features_lines

    alignment_progress$status_message <- "Traitement de la séquence ADN..."
    seq_lines <- gb_lines[(origin_line[1] + 1):length(gb_lines)]
    seq_raw <- paste(seq_lines, collapse = "")
    seq_clean <- gsub("[^acgtACGTnN]", "", seq_raw)
    data_xdna$seq <- Biostrings::DNAString(toupper(seq_clean))

    alignment_progress$status_message <- "Séquence chargée avec succès"
    Sys.sleep(0.5)
    alignment_progress$in_progress <- FALSE
    alignment_progress$complete <- TRUE

    showNotification("✅ Séquence GenBank chargée avec succès !",
                     type = "message", duration = 3)
  })

  seqs <- eventReactive(input$align_btn, {
    req(input$seq_files)

    lapply(input$seq_files, function(f) {
      lines <- readLines(f, warn = FALSE)
      seq_raw <- paste(lines, collapse = "")
      seq_clean <- clean_sequence(seq_raw)
      Biostrings::DNAString(seq_clean)
    })
  })

  # ==============================================================================
  # OUTPUTS PROGRESSION
  # ==============================================================================

  output$alignment_in_progress <- reactive({
    alignment_progress$in_progress
  })
  outputOptions(output, "alignment_in_progress", suspendWhenHidden = FALSE)

  output$processing_status <- renderUI({
    if (alignment_progress$in_progress) {
      tags$div(
        style = "color: #1976d2; font-weight: bold;",
        "⏳ ", alignment_progress$status_message
      )
    } else if (alignment_progress$complete) {
      tags$div(
        style = "color: #388e3c; font-weight: bold;",
        "✅ Prêt pour l'alignement"
      )
    } else {
      tags$div(
        style = "color: #757575;",
        "💤 En attente..."
      )
    }
  })

  output$alignment_complete <- reactive({
    if (alignment_progress$complete) {
      Sys.sleep(0.5)
      return(TRUE)
    }
    return(FALSE)
  })
  outputOptions(output, "alignment_complete", suspendWhenHidden = FALSE)

  # ==============================================================================
  # AFFICHAGE INFORMATIONS
  # ==============================================================================

  output$seq_xdna_compact <- renderText({
    req(data_xdna$seq)

    seq_text <- paste0("📋 Séquence GenBank (", length(data_xdna$seq), " nt):\n")
    seq_text <- paste0(seq_text, as.character(data_xdna$seq), "\n\n")

    seq_text <- paste0(seq_text, "🏷️ Annotations (Features):\n")
    if (!is.null(data_xdna$features) && length(data_xdna$features) > 0) {
      for (i in seq_along(data_xdna$features)) {
        seq_text <- paste0(seq_text, sprintf("%3d: %s\n", i, data_xdna$features[i]))
      }
    } else {
      seq_text <- paste0(seq_text, "Aucune annotation trouvée\n")
    }

    return(seq_text)
  })

  output$seqs_selected_compact <- renderText({
    req(seqs())

    selected_files <- input$seq_files
    if (!is.null(selected_files)) {
      result <- paste0("🧬 Séquences chargées (", length(seqs()), "):\n\n")

      for (i in seq_along(seqs())) {
        file_name <- basename(selected_files[i])
        seq_length <- length(seqs()[[i]])
        seq_string <- as.character(seqs()[[i]])

        result <- paste0(result, sprintf("=== %d. %s (%d nt) ===\n", i, file_name, seq_length))
        result <- paste0(result, seq_string, "\n\n")
      }
      return(result)
    } else {
      return("Aucune séquence sélectionnée")
    }
  })

  # ==============================================================================
  # GÉNÉRATION ALIGNEMENTS AVEC COULEURS
  # ==============================================================================

  output$align_results <- renderUI({
    req(data_xdna$seq, seqs())

    withProgress(message = 'Génération des alignements...', value = 0, {

      incProgress(0.1, detail = "Préparation des données...")

      # Génération de la légende des couleurs
      legend_content <- generate_color_legend(data_xdna$features, restriction_sites())

      align_output <- character()
      text_output <- character()
      selected_files <- input$seq_files

      total_seqs <- length(seqs())

      # Traitement de chaque séquence pour alignement
      for (i in seq_along(seqs())) {
        incProgress(0.7/total_seqs, detail = paste("Alignement", i, "sur", total_seqs))

        # Alignement par paires avec la séquence de référence
        aln <- Biostrings::pairwiseAlignment(
          pattern = seqs()[[i]],
          subject = data_xdna$seq,
          type = "local",
          substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1),
          gapOpening = -2,
          gapExtension = -1
        )

        # Extraction des informations d'alignement
        aln_start <- start(subject(aln))
        aln_end <- end(subject(aln))
        pat_aligned <- as.character(pattern(aln))
        sub_aligned <- as.character(subject(aln))
        annot_aligned <- annotate_sequence_mutations(pat_aligned, sub_aligned)

        # Création des séquences complètes avec gaps
        full_pattern <- create_full_pattern_with_gaps(pat_aligned, aln_start, length(data_xdna$seq))
        full_subject <- as.character(data_xdna$seq)
        full_annot <- create_full_annotation_with_spaces(annot_aligned, aln_start, length(data_xdna$seq))

        # Préparation des cartes de couleurs et restrictions
        colors <- build_sequence_color_map(data_xdna$features, 1, length(data_xdna$seq))
        restriction_positions <- build_restriction_position_map(length(data_xdna$seq), restriction_sites())

        # Nom d'affichage du fichier
        file_display_name <- if (!is.null(selected_files) && i <= length(selected_files)) {
          basename(selected_files[i])
        } else {
          paste0("Fichier_", i)
        }

        # Génération de l'alignement coloré
        blast_alignment <- generate_colored_alignment(
          full_pattern, full_subject, full_annot, aln_start, colors,
          paste0(file_display_name, " (région alignée: ", aln_start, "-", aln_end, ")"),
          restriction_positions
        )

        align_output <- c(align_output, blast_alignment$html)
        text_output <- c(text_output, blast_alignment$text)
      }

      incProgress(0.2, detail = "Finalisation...")

      # Sauvegarde pour les exports
      alignment_data$results <- align_output
      alignment_data$text_version <- paste(text_output, collapse = "")

      # Structure HTML finale avec légende et alignements côte à côte
      tags$div(
        id = "align_results",
        class = "results-container",
        tags$div(
          class = "legend-container",
          HTML(legend_content)
        ),
        tags$div(
          class = "alignments-container",
          HTML(paste(align_output, collapse = ""))
        )
      )
    })
  })

  # ==============================================================================
  # TÉLÉCHARGEMENTS COMPLETS
  # ==============================================================================

  output$download_html <- downloadHandler(
    filename = function() {
      paste0("alignement_complet_", Sys.Date(), "_", format(Sys.time(), "%H%M"), ".html")
    },
    content = function(file) {
      req(data_xdna$seq)

      selected_files <- input$seq_files
      files_display <- if (!is.null(selected_files)) {
        paste(basename(selected_files), collapse = ", ")
      } else {
        "Aucun fichier sélectionné"
      }

      html_content <- paste0(
        "<!DOCTYPE html>\n<html lang='fr'>\n<head>\n",
        "<meta charset='UTF-8'>\n",
        "<title>Rapport d'alignement HGX - ", Sys.Date(), "</title>\n",
        "<style>\n",
        "body { font-family: Arial, sans-serif; margin: 20px; background: #f8f9fa; }\n",
        ".container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; }\n",
        ".sequence-box { background: #ffffff; padding: 15px; border: 1px solid #dee2e6; border-radius: 4px; font-family: 'Courier New', monospace; font-size: 10px; white-space: pre-wrap; word-wrap: break-word; margin: 10px 0; }\n",
        ".legend-item { margin: 5px 0; }\n",
        ".color-box { display: inline-block; width: 16px; height: 16px; margin-right: 8px; border: 1px solid #ccc; }\n",
        "h1 { color: #b22222; margin: 0; }\n",
        "h2 { color: #2c3e50; border-bottom: 2px solid #dee2e6; padding-bottom: 10px; }\n",
        ".alignment-block { margin: 20px 0; padding: 15px; background: #ffffff; border: 1px solid #e0e0e0; border-radius: 5px; }\n",
        ".alignment-title { color: #2c3e50; font-weight: bold; margin-bottom: 10px; border-bottom: 2px solid #b22222; padding-bottom: 5px; }\n",
        "table { border-collapse: collapse; font-family: 'Courier New', monospace; }\n",
        "td { padding: 0; margin: 0; text-align: center; width: 1ch; min-width: 1ch; max-width: 1ch; }\n",
        "</style>\n</head>\n<body>\n"
      )

      html_content <- paste0(html_content, "<div class='container'>\n")
      html_content <- paste0(html_content, "<h1>🧬 RAPPORT D'ALIGNEMENT HGX</h1>\n")
      html_content <- paste0(html_content, "<p>Généré le ", Sys.time(), "</p>\n")
      html_content <- paste0(html_content, "<p><strong>Carte GenBank:</strong> ", input$carte_xdna, "</p>\n")
      html_content <- paste0(html_content, "<p><strong>Fichiers:</strong> ", files_display, "</p>\n")

      html_content <- paste0(html_content, "<h2>📋 Séquence GenBank de référence</h2>\n")
      html_content <- paste0(html_content, "<div class='sequence-box'>", as.character(data_xdna$seq), "</div>\n")

      if (!is.null(selected_files) && length(seqs()) > 0) {
        html_content <- paste0(html_content, "<h2>🧬 Séquences testées</h2>\n")
        for (i in seq_along(seqs())) {
          file_name <- basename(selected_files[i])
          seq_string <- as.character(seqs()[[i]])
          html_content <- paste0(html_content, "<h3>", i, ". ", file_name, "</h3>\n")
          html_content <- paste0(html_content, "<div class='sequence-box'>", seq_string, "</div>\n")
        }
      }

      # Légende des couleurs
      feats <- parse_genbank_features(data_xdna$features)
      if (length(feats) > 0) {
        html_content <- paste0(html_content, "<h2>🎨 Légende des couleurs</h2>\n")
        for (feat in feats) {
          if (!is.na(feat$position_raw)) {
            bounds <- as.numeric(unlist(strsplit(feat$position_raw, "\\.\\.")))
            if (length(bounds) == 2) {
              name_display <- if (feat$name != "") feat$name else paste("Feature", feat$type)
              color_to_use <- feat$color
              if (color_to_use == "#000000") {
                color_to_use <- get_color_by_feature_name(feat$name)
              }
              html_content <- paste0(html_content,
                                     "<div class='legend-item'>",
                                     "<span class='color-box' style='background-color:", color_to_use, ";'></span>",
                                     "<strong>", name_display, "</strong> (", bounds[1], "-", bounds[2], ")",
                                     "</div>\n")
            }
          }
        }
      }

      # Sites de restriction
      sites <- restriction_sites()
      if (length(sites) > 0) {
        html_content <- paste0(html_content, "<h3>Sites de restriction</h3>\n")
        enzymes <- get_restriction_enzymes()
        colors <- c("#FF0000", "#0000FF", "#00FF00", "#FF8000", "#8000FF", "#00FFFF")
        site_index <- 1

        for (enzyme_name in names(sites)) {
          enzyme_sites <- sites[[enzyme_name]]
          enzyme_seq <- enzymes[[enzyme_name]]
          color <- colors[((site_index - 1) %% length(colors)) + 1]

          html_content <- paste0(html_content,
                                 "<div class='legend-item'>",
                                 "<span class='color-box' style='background-color:", color, ";'></span>",
                                 "<strong>", enzyme_name, "</strong> (", enzyme_seq, ") - Sites: ", paste(enzyme_sites, collapse = ", "),
                                 "</div>\n")
          site_index <- site_index + 1
        }
      }

      # Alignements avec couleurs
      if (!is.null(alignment_data$results)) {
        clean_results <- alignment_data$results
        clean_results <- gsub('<span class="tooltip">[^<]*</span>', '', clean_results)
        clean_results <- gsub('class="nucleotide-cell"', '', clean_results)

        html_content <- paste0(html_content, "<h2>🔬 Alignements détaillés</h2>\n")
        html_content <- paste0(html_content, "<div style='overflow-x: auto;'>")
        html_content <- paste0(html_content, paste(clean_results, collapse = "\n"))
        html_content <- paste0(html_content, "</div>\n")
      }

      html_content <- paste0(html_content, "</div>\n</body>\n</html>")

      writeLines(html_content, file, useBytes = TRUE)
    }
  )

  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("alignement_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      req(data_xdna$seq, seqs())

      fasta_content <- character()

      fasta_content <- c(fasta_content,
                         paste0(">Sequence_Reference_", gsub("\\.gb$", "", input$carte_xdna)),
                         as.character(data_xdna$seq))

      selected_files <- input$seq_files
      for (i in seq_along(seqs())) {
        if (!is.null(selected_files) && i <= length(selected_files)) {
          seq_name <- gsub("\\.seq$", "", basename(selected_files[i]))
        } else {
          seq_name <- paste0("Unknown_", i)
        }
        fasta_content <- c(fasta_content,
                           paste0(">Sequence_", seq_name),
                           as.character(seqs()[[i]]))
      }

      writeLines(fasta_content, file)
    }
  )

}
