# ==============================================================================
# SERVER_CLONAGE.R
# Logique serveur pour l'application Shiny HGX - Module Clonage
# Gestion de la recherche de fichiers, alignement de séquences et analyse
# ==============================================================================

library(shiny)
library(Biostrings)

server_clonage <- function(input, output, session) {

  # ==============================================================================
  # VARIABLES RÉACTIVES
  # ==============================================================================

  # Données de la séquence de référence (GenBank)
  data_xdna <- reactiveValues(
    seq = NULL,        # Séquence ADN de référence
    features = NULL    # Annotations GenBank
  )

  # Données d'alignement
  alignment_data <- reactiveValues(
    results = NULL,      # Résultats HTML des alignements
    text_version = NULL  # Version texte pour export
  )

  # Données de recherche de fichiers
  search_data <- reactiveValues(
    folders_found = character(),      # Dossiers trouvés
    seq_files_found = character(),    # Noms d'affichage des fichiers .seq
    seq_files_paths = character(),    # Chemins complets des fichiers .seq
    search_in_progress = FALSE,       # Indicateur de recherche en cours
    stop_search = FALSE               # Signal d'arrêt de recherche
  )

  # ==============================================================================
  # CONFIGURATION DES CHEMINS
  # ==============================================================================

  # Répertoires de travail
  xdna_dir <- "P:/SEQ/Atest_cae"    # Dossier des fichiers GenBank
  seq_base_dir <- "P:/SEQ"          # Dossier racine des séquences

  # ==============================================================================
  # FONCTIONS DE RECHERCHE ULTRA-RAPIDE
  # ==============================================================================

  #' Recherche ultra-rapide de dossiers par pattern
  #' Utilise Sys.glob pour une recherche directe très performante
  #' @param plate_keyword Mot-clé de recherche pour la plaque
  #' @param base_dir Répertoire de base pour la recherche
  #' @return Chemin du dossier le plus récent trouvé ou vecteur vide
  search_ultra_fast <- function(plate_keyword, base_dir = "P:/SEQ") {
    if (is.null(plate_keyword) || plate_keyword == "") {
      return(character())
    }

    # Validation du répertoire de base
    if (!dir.exists(base_dir)) {
      return(character())
    }

    # Interface de progression Shiny
    withProgress(message = 'Recherche ultra-rapide...', value = 0, {

      tryCatch({
        # Nettoyage et préparation
        plate_keyword <- trimws(plate_keyword)
        start_time <- Sys.time()

        incProgress(0.3, detail = paste("Recherche directe de '*", plate_keyword, "*'...", sep=""))

        # Recherche directe par pattern avec Sys.glob (très rapide)
        search_pattern <- paste0("*", plate_keyword, "*")
        matching_folders <- Sys.glob(file.path(base_dir, search_pattern))

        # Filtrage pour ne garder que les dossiers
        matching_folders <- matching_folders[file.info(matching_folders)$isdir]
        search_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

        if (length(matching_folders) == 0) {
          incProgress(1, detail = paste("Aucun dossier trouvé en", round(search_time, 2), "s"))
          print(paste("❌ Aucun dossier trouvé pour le pattern:", search_pattern))
          return(character())
        }

        incProgress(0.3, detail = paste(length(matching_folders), "dossiers trouvés, tri..."))

        # Tri par nom décroissant pour avoir le plus récent en premier
        matching_folders <- sort(matching_folders, decreasing = TRUE)
        selected_folder <- matching_folders[1]

        total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        incProgress(0.4, detail = paste("Sélectionné:", basename(selected_folder), "en", round(total_time, 2), "s"))

        # Logs de diagnostic
        print(paste("✅ Recherche terminée en", round(total_time, 2), "secondes"))
        print(paste("📂 Dossiers trouvés:", length(matching_folders)))
        print(paste("🏆 Dossier sélectionné:", basename(selected_folder)))

        return(selected_folder)

      }, error = function(e) {
        incProgress(1, detail = paste("Erreur:", e$message))
        print(paste("❌ Erreur:", e$message))
        return(character())
      })
    })
  }

  #' Recherche exhaustive de fichiers .seq dans un dossier
  #' @param folder_path Chemin du dossier à analyser
  #' @param seq_keyword Mot-clé pour filtrer les fichiers .seq
  #' @return Liste avec noms d'affichage et chemins complets des fichiers
  search_all_seq_in_folder <- function(folder_path, seq_keyword) {
    if (length(folder_path) == 0 || is.null(seq_keyword) || seq_keyword == "") {
      return(list(files = character(), paths = character()))
    }

    tryCatch({
      seq_keyword <- trimws(seq_keyword)
      print(paste("🔍 Recherche fichiers .seq avec mot-clé:", seq_keyword))
      print(paste("📁 Dans le dossier:", basename(folder_path)))

      # Validation du dossier
      if (!dir.exists(folder_path)) {
        print(paste("❌ Dossier n'existe pas:", folder_path))
        return(list(files = character(), paths = character()))
      }

      # Recherche récursive de tous les fichiers .seq
      print("🔄 Scan de tous les fichiers .seq...")
      seq_files <- list.files(folder_path, pattern = "\\.seq$",
                              recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

      print(paste("📊 Total fichiers .seq trouvés:", length(seq_files)))

      if (length(seq_files) == 0) {
        print("❌ Aucun fichier .seq trouvé dans ce dossier")
        return(list(files = character(), paths = character()))
      }

      # Filtrage par mot-clé dans le nom de fichier
      file_names <- basename(seq_files)
      matching_indices <- grep(seq_keyword, file_names, ignore.case = TRUE)

      if (length(matching_indices) == 0) {
        print(paste("❌ Aucun fichier .seq correspondant au mot-clé '", seq_keyword, "'"))
        # Affichage d'exemples pour aider l'utilisateur
        print(paste("📋 Exemples de fichiers trouvés:",
                    paste(file_names[1:min(5, length(file_names))], collapse = ", ")))
        return(list(files = character(), paths = character()))
      }

      print(paste("✅ Fichiers correspondants au mot-clé:", length(matching_indices)))

      matching_files <- seq_files[matching_indices]
      matching_names <- file_names[matching_indices]

      # Création de noms d'affichage avec chemin relatif
      folder_name <- basename(folder_path)
      display_names <- character()

      for (i in seq_along(matching_files)) {
        # Calcul du chemin relatif
        relative_path <- gsub(paste0("^", gsub("([\\(\\)\\[\\]\\{\\}\\^\\$\\*\\+\\?\\|\\\\])",
                                               "\\\\\\1", folder_path), "[\\/\\\\]?"), "", matching_files[i])

        # Formatage du nom d'affichage avec indication du sous-dossier si applicable
        if (dirname(relative_path) != ".") {
          display_names[i] <- paste0("[", folder_name, "] ", dirname(relative_path), " → ", basename(relative_path))
        } else {
          display_names[i] <- paste0("[", folder_name, "] ", basename(relative_path))
        }
      }

      print(paste("📋 Exemples de fichiers trouvés:",
                  paste(display_names[1:min(3, length(display_names))], collapse = ", ")))

      return(list(files = display_names, paths = matching_files))

    }, error = function(e) {
      print(paste("❌ Erreur dans search_all_seq_in_folder:", e$message))
      return(list(files = character(), paths = character()))
    })
  }

  # ==============================================================================
  # GESTION DES FICHIERS GENBANK
  # ==============================================================================

  #' Récupération de la liste des fichiers GenBank disponibles
  #' @return Vecteur des noms de fichiers .gb
  get_available_gb_files <- function() {
    gb_files <- list.files(xdna_dir, pattern = "\\.gb$", full.names = FALSE)
    return(gb_files)
  }

  # Initialisation de la liste des fichiers GenBank au démarrage
  observe({
    gb_files <- get_available_gb_files()
    updateSelectInput(session, "carte_xdna", choices = gb_files)
  })

  # Rafraîchissement manuel de la liste des fichiers GenBank
  observeEvent(input$refresh_files, {
    gb_files <- get_available_gb_files()
    updateSelectInput(session, "carte_xdna", choices = gb_files)
    showNotification("📁 Liste des fichiers GenBank mise à jour !", type = "message", duration = 2)
  })

  # ==============================================================================
  # ANALYSE DES SITES DE RESTRICTION
  # ==============================================================================

  #' Calcul réactif des sites de restriction
  #' Analyse la séquence de référence avec les enzymes sélectionnées
  restriction_sites <- reactive({
    req(data_xdna$seq)
    sites_list <- list()
    enzymes <- get_restriction_enzymes()

    # Analyse enzyme 1
    if (!is.null(input$enzyme1) && input$enzyme1 != "") {
      sites1 <- find_restriction_sites(data_xdna$seq, enzymes[[input$enzyme1]])
      if (length(sites1) > 0) sites_list[[input$enzyme1]] <- sites1
    }

    # Analyse enzyme 2
    if (!is.null(input$enzyme2) && input$enzyme2 != "") {
      sites2 <- find_restriction_sites(data_xdna$seq, enzymes[[input$enzyme2]])
      if (length(sites2) > 0) sites_list[[input$enzyme2]] <- sites2
    }

    return(sites_list)
  })

  # Affichage des informations sur les sites de restriction trouvés
  output$restriction_info <- renderText({
    sites <- restriction_sites()
    if (length(sites) == 0) return("Aucun site trouvé")

    info <- sapply(names(sites), function(e) paste0(e, ": ", length(sites[[e]]), " site(s)"))
    paste(info, collapse = " | ")
  })

  # ==============================================================================
  # GESTION DE LA RECHERCHE DE FICHIERS
  # ==============================================================================

  # Signal d'arrêt de recherche
  observeEvent(input$stop_search, {
    search_data$stop_search <- TRUE
    search_data$search_in_progress <- FALSE
    showNotification("🛑 Recherche interrompue", type = "warning", duration = 2)
  })

  # Événement principal de recherche de fichiers .seq
  observeEvent(input$search_seq_btn, {
    req(input$plate_keyword, input$seq_keyword)

    # Validation des entrées utilisateur
    if (nchar(trimws(input$plate_keyword)) == 0) {
      showNotification("⚠️ Veuillez saisir un nom de plaque", type = "warning", duration = 3)
      return()
    }

    if (nchar(trimws(input$seq_keyword)) == 0) {
      showNotification("⚠️ Veuillez saisir un mot-clé pour les fichiers .seq", type = "warning", duration = 3)
      return()
    }

    # Réinitialisation des flags et données
    search_data$stop_search <- FALSE
    search_data$search_in_progress <- TRUE
    search_data$folders_found <- character()
    search_data$seq_files_found <- character()
    search_data$seq_files_paths <- character()
    updateSelectInput(session, "seq_files", choices = NULL)

    # Processus de recherche principal
    tryCatch({
      # Phase 1: Recherche ultra-rapide du dossier de plaque
      plate_folder <- search_ultra_fast(input$plate_keyword, seq_base_dir)

      if (length(plate_folder) == 0) {
        search_data$search_in_progress <- FALSE
        showNotification("❌ Aucun dossier trouvé pour cette plaque", type = "error", duration = 5)
        return()
      }

      search_data$folders_found <- plate_folder
      showNotification(paste("📂 Dossier sélectionné:", basename(plate_folder)), type = "message", duration = 3)

      # Phase 2: Recherche exhaustive des fichiers .seq
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
            # Mise à jour de l'interface avec les fichiers trouvés
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
      print(paste("❌ Erreur de recherche:", e$message))
    })
  })

  # ==============================================================================
  # OUTPUTS POUR L'INTERFACE UTILISATEUR
  # ==============================================================================

  # Indicateur de recherche en cours (pour affichage conditionnel)
  output$search_in_progress <- reactive({
    search_data$search_in_progress
  })
  outputOptions(output, "search_in_progress", suspendWhenHidden = FALSE)

  # Affichage formaté des résultats de recherche
  output$search_results <- renderText({
    req(length(search_data$folders_found) > 0 || length(search_data$seq_files_found) > 0)

    if (length(search_data$folders_found) > 0) {
      folder_name <- basename(search_data$folders_found[1])
      folders_text <- paste0("📂 Dossier sélectionné: ", folder_name)

      if (length(search_data$seq_files_found) > 0) {
        # Limitation de l'affichage pour éviter la surcharge
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

  # Condition d'affichage du sélecteur de fichiers
  output$seq_files_found <- reactive({
    length(search_data$seq_files_found) > 0
  })
  outputOptions(output, "seq_files_found", suspendWhenHidden = FALSE)

  # ==============================================================================
  # CHARGEMENT ET TRAITEMENT DES FICHIERS GENBANK
  # ==============================================================================

  # Événement de chargement et alignement des séquences
  observeEvent(input$align_btn, {
    req(input$carte_xdna)
    fichier <- file.path(xdna_dir, input$carte_xdna)

    # Lecture robuste du fichier GenBank avec gestion d'encodage
    tryCatch({
      # Tentative de lecture en UTF-8
      gb_lines <- readLines(fichier, warn = FALSE, encoding = "UTF-8")
    }, error = function(e) {
      # Fallback sur latin1 si UTF-8 échoue
      tryCatch({
        gb_lines <- readLines(fichier, warn = FALSE, encoding = "latin1")
      }, error = function(e2) {
        # Dernière option: lecture brute
        gb_lines <- readLines(fichier, warn = FALSE)
      })
    })

    # Nettoyage des caractères problématiques
    gb_lines <- iconv(gb_lines, to = "UTF-8", sub = "")
    gb_lines <- gb_lines[!is.na(gb_lines)]
    gb_lines <- gsub("[^\x01-\x7F]", "", gb_lines)  # Suppression caractères non-ASCII

    # Localisation de la section ORIGIN (début de la séquence)
    origin_line <- tryCatch({
      grep("^ORIGIN", gb_lines, ignore.case = TRUE)
    }, warning = function(w) {
      which(grepl("^ORIGIN", gb_lines, ignore.case = TRUE))
    })

    if (length(origin_line) == 0) {
      showNotification("Erreur: Section ORIGIN non trouvée dans le fichier GenBank",
                       type = "error", duration = 5)
      return()
    }

    # Extraction des features (annotations)
    features_block <- gb_lines[1:(origin_line[1] - 1)]
    features_lines <- tryCatch({
      features_block[grep("^\\s{5}|^\\s{21}", features_block)]
    }, warning = function(w) {
      features_block[grepl("^\\s{5}|^\\s{21}", features_block)]
    })

    data_xdna$features <- features_lines

    # Extraction et nettoyage de la séquence ADN
    seq_lines <- gb_lines[(origin_line[1] + 1):length(gb_lines)]
    seq_raw <- paste(seq_lines, collapse = "")
    seq_clean <- gsub("[^acgtACGTnN]", "", seq_raw)
    data_xdna$seq <- Biostrings::DNAString(toupper(seq_clean))
  })

  # Chargement réactif des séquences .seq sélectionnées
  seqs <- eventReactive(input$align_btn, {
    req(input$seq_files)

    # Traitement de chaque fichier .seq sélectionné
    lapply(input$seq_files, function(f) {
      lines <- readLines(f, warn = FALSE)
      seq_raw <- paste(lines, collapse = "")
      seq_clean <- clean_sequence(seq_raw)
      Biostrings::DNAString(seq_clean)
    })
  })

  # ==============================================================================
  # AFFICHAGE DES INFORMATIONS
  # ==============================================================================

  # Affichage de la séquence de référence et de ses annotations
  output$seq_xdna <- renderPrint({
    req(data_xdna$seq)
    cat("Séquence GenBank :\n")
    print(data_xdna$seq)
    cat("\nAnnotations (Features) :\n")
    print(data_xdna$features)
  })

  # Affichage des séquences sélectionnées
  output$seqs_selected <- renderPrint({
    req(seqs())
    cat("Séquences chargées des fichiers .seq :\n")

    selected_files <- input$seq_files
    if (!is.null(selected_files)) {
      for (i in seq_along(seqs())) {
        file_name <- basename(selected_files[i])
        cat(" - ", file_name, ": ", as.character(seqs()[[i]]), "\n")
      }
    }
  })

  # ==============================================================================
  # GÉNÉRATION DES ALIGNEMENTS
  # ==============================================================================

  # Génération et affichage des résultats d'alignement
  output$align_results <- renderUI({
    req(data_xdna$seq, seqs())

    # Génération de la légende des couleurs
    legend_content <- generate_color_legend(data_xdna$features, restriction_sites())

    align_output <- character()
    text_output <- character()
    selected_files <- input$seq_files

    # Traitement de chaque séquence pour alignement
    for (i in seq_along(seqs())) {
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

  # ==============================================================================
  # GESTIONNAIRES DE TÉLÉCHARGEMENT
  # ==============================================================================

  # Téléchargement au format TXT avec en-têtes et légendes
  output$download_txt <- downloadHandler(
    filename = function() {
      paste0("alignement_", Sys.Date(), ".txt")
    },
    content = function(file) {
      req(alignment_data$text_version)

      # Préparation de la liste des fichiers pour l'en-tête
      selected_files <- input$seq_files
      files_display <- if (!is.null(selected_files)) {
        paste(basename(selected_files), collapse = ", ")
      } else {
        "Aucun fichier sélectionné"
      }

      # En-tête du fichier
      header <- paste0(
        "=== RESULTATS D'ALIGNEMENT - HGX ===\n",
        "Date: ", Sys.time(), "\n",
        "Carte GenBank: ", input$carte_xdna, "\n",
        "Fichiers séquences: ", files_display, "\n",
        "Longueur séquence référence: ", length(data_xdna$seq), " nt\n\n",
        "LEGENDE DES COULEURS:\n"
      )

      # Génération de la légende des features
      feats <- parse_genbank_features(data_xdna$features)
      legend_text <- ""
      for (feat in feats) {
        if (!is.na(feat$position_raw)) {
          bounds <- as.numeric(unlist(strsplit(feat$position_raw, "\\.\\.")))
          if (length(bounds) == 2) {
            name_display <- if (feat$name != "") feat$name else "Feature sans nom"
            legend_text <- paste0(legend_text, "- ", name_display, " (", bounds[1], "-", bounds[2],
                                  ") - Couleur: ", feat$color, "\n")
          }
        }
      }

      # Ajout de la légende des sites de restriction
      sites <- restriction_sites()
      if (length(sites) > 0) {
        legend_text <- paste0(legend_text, "\nSITES DE RESTRICTION:\n")
        enzymes <- get_restriction_enzymes()
        for (enzyme_name in names(sites)) {
          enzyme_sites <- sites[[enzyme_name]]
          enzyme_seq <- enzymes[[enzyme_name]]
          legend_text <- paste0(legend_text, "- ", enzyme_name, " (", enzyme_seq,
                                ") - Sites: ", paste(enzyme_sites, collapse = ", "), "\n")
        }
      }

      # Assemblage final et écriture
      content <- paste0(header, legend_text, "\n", alignment_data$text_version)
      writeLines(content, file, useBytes = TRUE)
    }
  )

  # Téléchargement au format FASTA
  output$download_fasta <- downloadHandler(
    filename = function() {
      paste0("alignement_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      req(data_xdna$seq, seqs())

      fasta_content <- character()

      # Ajout de la séquence de référence
      fasta_content <- c(fasta_content,
                         paste0(">Sequence_Reference_", gsub("\\.gb$", "", input$carte_xdna)),
                         as.character(data_xdna$seq))

      # Ajout des séquences testées
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
