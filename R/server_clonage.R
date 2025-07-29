# ==============================================================================
# SERVER_CLONAGE.R
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
    stop_search = FALSE,
    file_groups = list()
  )

  alignment_progress <- reactiveValues(
    in_progress = FALSE,
    complete = FALSE,
    status_message = "Prêt"
  )

  # ==============================================================================
  # CONFIGURATION DES CHEMINS
  # ==============================================================================

  # Utiliser la configuration centralisée avec gestion d'erreur
  xdna_dir <- NULL
  seq_base_dir <- NULL

  tryCatch({
    if (exists("get_config") && is.function(get_config)) {
      cat("🔧 Appel de get_config()...\n")
      config <- get_config()

      cat("📋 Résultat de get_config():\n")
      cat("   - config$xdna_dir:", deparse(config$xdna_dir), "\n")
      cat("   - config$seq_dir:", deparse(config$seq_dir), "\n")
      cat("   - config$environment:", deparse(config$environment), "\n")

      # CORRECTION : Assignation directe sans <<-
      if (!is.null(config$xdna_dir) && config$xdna_dir != "") {
        xdna_dir <- config$xdna_dir
      }
      if (!is.null(config$seq_dir) && config$seq_dir != "") {
        seq_base_dir <- config$seq_dir
      }

      cat("📁 Variables assignées:\n")
      cat("   - xdna_dir:", xdna_dir, "\n")
      cat("   - seq_base_dir:", seq_base_dir, "\n")

    } else {
      stop("get_config non disponible")
    }

  }, error = function(e) {
    cat("⚠️ Erreur avec get_config():", e$message, "\n")
    cat("   Utilisation du fallback\n")

    # Configuration directe (fallback)
    if (.Platform$OS.type == "windows") {
      if (dir.exists("R:/Production/Labo YEAST/Demandes du service/carte_nouveaux_clonages")) {
        xdna_dir <- "R:/Production/Labo YEAST/Demandes du service/carte_nouveaux_clonages"
        seq_base_dir <- "P:/SEQ"
        cat("   - Mode: Windows développement\n")
      }
    } else {
      if (dir.exists("/mnt/carte_nouveaux_clonages")) {
        xdna_dir <- "/mnt/carte_nouveaux_clonages"
        seq_base_dir <- "/data/production/SEQ"
        cat("   - Mode: Linux production\n")
      }
    }
  })

  # Vérification finale et fallback ultime
  if (is.null(xdna_dir) || is.null(seq_base_dir) || xdna_dir == "" || seq_base_dir == "") {
    cat("❌ Variables toujours NULL, fallback ultime\n")

    if (.Platform$OS.type == "windows") {
      xdna_dir <- "R:/Production/Labo YEAST/Demandes du service/carte_nouveaux_clonages"
      seq_base_dir <- "P:/SEQ"
      cat("   - Fallback Windows forcé\n")
    } else {
      xdna_dir <- "/mnt/carte_nouveaux_clonages"
      seq_base_dir <- "/data/production/SEQ"
      cat("   - Fallback Linux forcé\n")
    }
  }

  cat("✅ Configuration finale:\n")
  cat("   - xdna_dir:", xdna_dir, "\n")
  cat("   - seq_base_dir:", seq_base_dir, "\n")

  # Test d'accès aux fichiers
  cat("\n🔍 TEST D'ACCÈS AUX FICHIERS:\n")
  cat("Dossier testé:", xdna_dir, "\n")

  if (dir.exists(xdna_dir)) {
    cat("✅ Dossier accessible\n")

    # Lister tous les fichiers
    all_files <- list.files(xdna_dir, full.names = FALSE)
    cat("   - Total fichiers:", length(all_files), "\n")

    if (length(all_files) > 0) {
      cat("   - Premiers fichiers:", paste(head(all_files, 3), collapse = ", "), "\n")

      # Chercher les .gb
      gb_files <- all_files[grepl("\\.gb$", all_files, ignore.case = TRUE)]
      cat("   - Fichiers .gb trouvés:", length(gb_files), "\n")

      if (length(gb_files) > 0) {
        cat("   - Fichiers .gb:", paste(gb_files, collapse = ", "), "\n")
      }
    }
  } else {
    cat("❌ Dossier non accessible:", xdna_dir, "\n")
  }

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
      return(list(files = character(), paths = character(), groups = list()))
    }

    tryCatch({
      seq_keyword <- trimws(seq_keyword)

      if (!dir.exists(folder_path)) {
        return(list(files = character(), paths = character(), groups = list()))
      }

      seq_files <- list.files(folder_path, pattern = "\\.seq$",
                              recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

      if (length(seq_files) == 0) {
        return(list(files = character(), paths = character(), groups = list()))
      }

      file_names <- basename(seq_files)
      matching_indices <- grep(seq_keyword, file_names, ignore.case = TRUE)

      if (length(matching_indices) == 0) {
        return(list(files = character(), paths = character(), groups = list()))
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

      # CHANGEMENT : Organiser par clones d'abord
      clones <- organize_files_by_fragments(matching_files, display_names)

      return(list(files = display_names, paths = matching_files, groups = clones))

    }, error = function(e) {
      return(list(files = character(), paths = character(), groups = list()))
    })
  }

  # ==============================================================================
  # GESTION DES FICHIERS GENBANK
  # ==============================================================================

  get_available_gb_files <- function() {
    if (is.null(xdna_dir) || !dir.exists(xdna_dir)) {
      cat("⚠️ get_available_gb_files: Dossier non accessible\n")
      return(character())
    }

    tryCatch({
      all_files <- list.files(xdna_dir, full.names = FALSE)
      gb_files <- all_files[grepl("\\.gb$", all_files, ignore.case = TRUE)]
      cat("📄 get_available_gb_files: ", length(gb_files), " fichiers .gb trouvés\n")
      return(gb_files)
    }, error = function(e) {
      cat("❌ Erreur get_available_gb_files:", e$message, "\n")
      return(character())
    })
  }

  # Initialisation sécurisée
  observe({
    tryCatch({
      gb_files <- get_available_gb_files()
      if (length(gb_files) == 0) {
        gb_files <- c("Aucun fichier .gb trouvé" = "")
        showNotification("⚠️ Aucun fichier GenBank trouvé. Vérifiez les montages.",
                         type = "warning", duration = 10)
      }
      updateSelectInput(session, "carte_xdna", choices = gb_files)
    }, error = function(e) {
      cat("❌ Erreur lors de l'initialisation des fichiers GB:", e$message, "\n")
      updateSelectInput(session, "carte_xdna", choices = c("Erreur de chargement" = ""))
      showNotification("❌ Erreur lors du chargement des fichiers GenBank",
                       type = "error", duration = 10)
    })
  })
  # ==============================================================================
  # SITES DE RESTRICTION
  # ==============================================================================

  restriction_sites <- reactive({
    req(data_xdna$seq)

    # Récupérer les informations sur les fragments si disponibles
    selected_files <- get_all_selected_files()
    if (length(selected_files) == 0) {
      return(list())
    }

    sites_list <- list()
    enzymes <- get_restriction_enzymes()

    # Déterminer les séquences d'enzymes
    enzyme1_seq <- NULL
    enzyme2_seq <- NULL

    if (!is.null(input$enzyme1) && input$enzyme1 != "") {
      if (input$enzyme1 == "CUSTOM") {
        enzyme1_seq <- validate_enzyme_sequence(input$enzyme1_custom_seq)
      } else {
        enzyme1_seq <- enzymes[[input$enzyme1]]
      }
    }

    if (!is.null(input$enzyme2) && input$enzyme2 != "") {
      if (input$enzyme2 == "CUSTOM") {
        enzyme2_seq <- validate_enzyme_sequence(input$enzyme2_custom_seq)
      } else {
        enzyme2_seq <- enzymes[[input$enzyme2]]
      }
    }

    # Chercher les enzymes dans la séquence de référence
    if (!is.null(enzyme1_seq)) {
      enzyme1_name <- if (input$enzyme1 == "CUSTOM") {
        if (!is.null(input$enzyme1_custom_name) && input$enzyme1_custom_name != "") {
          input$enzyme1_custom_name
        } else {
          paste0("Custom1_", enzyme1_seq)
        }
      } else {
        input$enzyme1
      }

      sites1 <- find_restriction_sites(data_xdna$seq, enzyme1_seq)
      if (length(sites1) > 0) {
        sites_list[[enzyme1_name]] <- sites1
        attr(sites_list[[enzyme1_name]], "enzyme_sequence") <- enzyme1_seq
      }
    }

    if (!is.null(enzyme2_seq)) {
      enzyme2_name <- if (input$enzyme2 == "CUSTOM") {
        if (!is.null(input$enzyme2_custom_name) && input$enzyme2_custom_name != "") {
          input$enzyme2_custom_name
        } else {
          paste0("Custom2_", enzyme2_seq)
        }
      } else {
        input$enzyme2
      }

      sites2 <- find_restriction_sites(data_xdna$seq, enzyme2_seq)
      if (length(sites2) > 0) {
        sites_list[[enzyme2_name]] <- sites2
        attr(sites_list[[enzyme2_name]], "enzyme_sequence") <- enzyme2_seq
      }
    }

    return(sites_list)
  })


  output$restriction_info <- renderText({
    sites <- restriction_sites()
    selected_files <- get_all_selected_files()

    if (length(sites) == 0) {
      return("Aucun site trouvé dans la séquence de référence")
    }

    info_parts <- character()

    # Informations sur les sites trouvés
    for (enzyme_name in names(sites)) {
      enzyme_sites <- sites[[enzyme_name]]
      info_parts <- c(info_parts, paste0(enzyme_name, ": ", length(enzyme_sites), " site(s)"))
    }

    # Informations sur les fragments si disponibles
    if (length(selected_files) > 0) {
      fragment_types <- sapply(selected_files, extract_fragment_type)
      fragment_summary <- table(fragment_types)

      if (length(fragment_summary) > 0) {
        fragment_info <- paste(names(fragment_summary), fragment_summary, sep = ":", collapse = " | ")
        info_parts <- c(info_parts, paste0("Fragments: ", fragment_info))
      }
    }

    return(paste(info_parts, collapse = " | "))
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
    search_data$file_groups <- list()  # Réinitialiser les groupes


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
          search_data$file_groups <- seq_results$groups  # Stocker les groupes

          incProgress(0.5, detail = paste(length(seq_results$files), "fichiers trouvés"))
          search_data$search_in_progress <- FALSE

          if (length(seq_results$files) == 0) {
            showNotification("❌ Aucun fichier .seq trouvé avec ce mot-clé", type = "warning", duration = 5)
          } else {
            showNotification(paste("✅", length(seq_results$files), "fichier(s) .seq trouvé(s) dans", length(seq_results$groups), "groupe(s)"),
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

  # Observer 1 : Pour les checkboxes individuelles des clones
  observe({
    req(search_data$file_groups)

    clones <- search_data$file_groups

    # Vérification que clones n'est pas vide
    if (length(clones) == 0) {
      return()
    }

    global_counter <- 0

    for (clone_id in names(clones)) {
      global_counter <- global_counter + 1
      clone_data <- clones[[clone_id]]

      # Vérification que clone_data existe et a des paths
      if (is.null(clone_data) || is.null(clone_data$paths)) {
        next
      }

      local({
        counter <- global_counter
        paths <- clone_data$paths

        observeEvent(input[[paste0("select_clone_", counter)]], {
          checkbox_value <- input[[paste0("select_clone_", counter)]]
          clone_select_id <- paste0("seq_files_clone_", counter)

          if (isTRUE(checkbox_value)) {
            updateSelectInput(session, clone_select_id, selected = paths)
          } else {
            updateSelectInput(session, clone_select_id, selected = character())
          }
        }, ignoreInit = TRUE, ignoreNULL = TRUE)
      })
    }
  })

  # ==============================================================================
  # AJOUT DANS SERVER_CLONAGE.R - VALIDATION EN TEMPS RÉEL
  # ==============================================================================

  # Validation séquence 1 personnalisée
  observeEvent(input$enzyme1_custom_seq, {
    if (!is.null(input$enzyme1_custom_seq) && input$enzyme1_custom_seq != "") {
      validated <- validate_enzyme_sequence(input$enzyme1_custom_seq)
      if (is.null(validated)) {
        showNotification("⚠️ Séquence 1 invalide. Utilisez seulement A, T, C, G (minimum 3 caractères)",
                         type = "warning", duration = 3)
      } else {
        showNotification("✅ Séquence 1 valide - Recherche possible",
                         type = "message", duration = 2)
      }
    }
  })

  # Validation séquence 2 personnalisée
  observeEvent(input$enzyme2_custom_seq, {
    if (!is.null(input$enzyme2_custom_seq) && input$enzyme2_custom_seq != "") {
      validated <- validate_enzyme_sequence(input$enzyme2_custom_seq)
      if (is.null(validated)) {
        showNotification("⚠️ Séquence 2 invalide. Utilisez seulement A, T, C, G (minimum 3 caractères)",
                         type = "warning", duration = 3)
      } else {
        showNotification("✅ Séquence 2 valide - Recherche possible",
                         type = "message", duration = 2)
      }
    }
  })

  # ==============================================================================
  # GESTION DES FICHIERS AB1
  # ==============================================================================

  # Variable réactive pour stocker les fichiers AB1
  ab1_data <- reactiveValues(
    file_info = list(),
    scanned = FALSE
  )

  # Scan automatique des fichiers AB1 quand la sélection change
  observe({
    selected_files <- get_all_selected_files()

    if (length(selected_files) > 0) {
      # Utiliser la fonction du global_clonage.R
      ab1_info <- find_corresponding_ab1_files(selected_files)

      # Stocker les informations
      ab1_data$file_info <- ab1_info
      ab1_data$scanned <- TRUE
    } else {
      # Réinitialiser si aucun fichier sélectionné
      ab1_data$file_info <- list()
      ab1_data$scanned <- FALSE
    }
  })

  # Observer unique pour AB1 - VERSION CORRIGÉE
  observe({
    req(ab1_data$file_info)
    is_server <- !(.Platform$OS.type == "windows")

    if (is_server) {
      # SERVEUR : Créer les downloadHandlers simples
      lapply(seq_along(ab1_data$file_info), function(i) {
        file_info <- ab1_data$file_info[[i]]

        if (file_info$exists) {
          output[[paste0("download_ab1_", i)]] <- downloadHandler(
            filename = function() {
              basename(file_info$ab1_file)
            },
            content = function(file) {
              file.copy(file_info$ab1_file, file)
            },
            contentType = "application/x-abi"
          )
        }
      })
    } else {
      # LOCAL : Observers pour ouverture directe
      lapply(seq_along(ab1_data$file_info), function(i) {
        file_info <- ab1_data$file_info[[i]]

        if (file_info$exists) {
          observeEvent(input[[paste0("open_ab1_", i)]], {
            file_path <- file_info$ab1_file
            file_name <- basename(file_path)

            success <- open_file_with_default_app(file_path)

            if (success) {
              showNotification(
                paste0("📂 Ouverture de ", file_name, "..."),
                type = "message", duration = 3
              )
            } else {
              showNotification(
                paste0("❌ Impossible d'ouvrir ", file_name),
                type = "error", duration = 5
              )
            }
          }, ignoreInit = TRUE, autoDestroy = TRUE)
        }
      })
    }
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
    !is.null(search_data$file_groups) && length(search_data$file_groups) > 0
  })
  outputOptions(output, "seq_files_found", suspendWhenHidden = FALSE)

  # Output pour l'interface des groupes
  output$groups_selection_ui <- renderUI({
    req(search_data$file_groups)

    clones <- search_data$file_groups

    if (length(clones) == 0) {
      return(NULL)
    }

    clone_uis <- list()

    # Compteur global pour les IDs uniques
    global_counter <- 0

    for (clone_id in names(clones)) {
      global_counter <- global_counter + 1
      clone_data <- clones[[clone_id]]

      # Vérifications de sécurité
      if (is.null(clone_data) || is.null(clone_data$paths)) {
        next
      }

      fragments <- clone_data$fragments

      # Calculer le nombre total de fichiers dans ce clone
      total_files <- length(clone_data$paths)

      # CORRECTION 1 : Créer la liste des choix TRIÉE dans l'ordre correct
      # Récupérer les paths et display_names dans l'ordre des fragments
      ordered_paths <- character()
      ordered_display_names <- character()

      # Ordre souhaité : 5p → int → 3p
      fragment_order <- c("5p", "int", "int1", "int2", "int3", "int4", "int5", "3p")

      for (fragment_type in fragment_order) {
        if (fragment_type %in% names(fragments)) {
          fragment_data <- fragments[[fragment_type]]
          if (!is.null(fragment_data) && !is.null(fragment_data$paths)) {
            ordered_paths <- c(ordered_paths, fragment_data$paths)
            ordered_display_names <- c(ordered_display_names, fragment_data$display_names)
          }
        }
      }

      # Si pas de tri possible, utiliser l'ordre original
      if (length(ordered_paths) == 0) {
        ordered_paths <- clone_data$paths
        ordered_display_names <- clone_data$display_names
      }

      choices_list <- setNames(ordered_paths, ordered_display_names)

      # Créer un résumé des fragments présents
      fragment_summary <- character()

      if (length(fragments) > 0) {
        for (fragment_type in names(fragments)) {
          fragment_data <- fragments[[fragment_type]]

          if (!is.null(fragment_data) && !is.null(fragment_data$paths)) {
            fragment_count <- length(fragment_data$paths)

            # Icône selon le type
            if (fragment_type == "5p") {
              icon <- "🔹"
            } else if (fragment_type == "3p") {
              icon <- "🔸"
            } else if (grepl("int", fragment_type)) {
              icon <- "🔻"
            } else {
              icon <- "🧬"
            }

            fragment_summary <- c(fragment_summary, paste0(icon, fragment_type, ":", fragment_count))
          }
        }
      }

      # Si pas de résumé, en créer un par défaut
      if (length(fragment_summary) == 0) {
        fragment_summary <- paste0("📁 ", total_files, " fichier(s)")
      }

      # Créer l'UI pour ce clone
      clone_ui <- div(
        style = "margin-bottom: 15px; padding: 12px; border: 2px solid #6c757d; border-radius: 6px; background: #f8f9fa;",

        # En-tête du clone
        div(style = "margin-bottom: 10px; padding-bottom: 8px; border-bottom: 1px solid #dee2e6;",
            div(style = "display: flex; justify-content: space-between; align-items: center;",
                div(
                  h5(paste("🧬", clone_id, "(", total_files, "fichiers)"),
                     style = "color: #495057; margin: 0; font-weight: bold;"),
                  tags$small(paste(fragment_summary, collapse = " | "),
                             style = "color: #6c757d; font-style: italic;")
                ),
                div(
                  checkboxInput(paste0("select_clone_", global_counter),
                                label = "Tout sélectionner",
                                value = FALSE)
                )
            )
        ),

        # Liste de sélection pour ce clone
        selectInput(paste0("seq_files_clone_", global_counter),
                    label = NULL,
                    choices = choices_list,
                    multiple = TRUE,
                    width = "100%",
                    selectize = FALSE,
                    size = min(6, length(choices_list)))
      )

      clone_uis[[clone_id]] <- clone_ui
    }

    # CORRECTION 2 : Supprimer le message vert
    return(div(
      h5("🧬 Clones trouvés", style = "color: #b22222; margin-bottom: 15px;"),
      clone_uis
    ))
  })

  # Fonction pour récupérer tous les fichiers sélectionnés
  get_all_selected_files <- reactive({
    req(search_data$file_groups)

    all_selected <- character()
    clones <- search_data$file_groups

    if (length(clones) == 0) {
      return(all_selected)
    }

    global_counter <- 0

    for (clone_id in names(clones)) {
      global_counter <- global_counter + 1
      clone_input_name <- paste0("seq_files_clone_", global_counter)
      clone_selection <- input[[clone_input_name]]

      if (!is.null(clone_selection) && length(clone_selection) > 0) {
        all_selected <- c(all_selected, clone_selection)
      }
    }

    return(all_selected)
  })


  # Output pour le résumé de sélection
  output$selection_summary_text <- renderText({
    selected_files <- get_all_selected_files()
    if (length(selected_files) == 0) {
      return("Aucun fichier sélectionné")
    } else {
      return(paste("📊 Total sélectionné:", length(selected_files), "fichier(s)"))
    }
  })

  # Interface dynamique pour les boutons AB1 individuels
  output$ab1_buttons_ui <- renderUI({
    req(ab1_data$scanned)

    if (length(ab1_data$file_info) == 0) {
      return(div(
        style = "padding: 20px; text-align: center; color: #6c757d;",
        "📝 Aucun fichier AB1 trouvé."
      ))
    }

    button_list <- list()
    is_server <- !(.Platform$OS.type == "windows")

    for (i in seq_along(ab1_data$file_info)) {
      file_info <- ab1_data$file_info[[i]]
      ab1_name <- basename(file_info$ab1_file)

      if (file_info$exists) {
        if (is_server) {
          # SERVEUR : Bouton de téléchargement simple
          button_list[[i]] <- div(
            style = "margin-bottom: 8px; padding: 8px; background: #f8f9fa; border: 1px solid #28a745; border-radius: 4px;",

            div(style = "display: flex; justify-content: space-between; align-items: center;",
                div(style = "flex: 1;",
                    tags$strong(style = "color: #28a745;", "✅ ", ab1_name),
                    br(),
                    tags$small(style = "color: #6c757d;", "Clic → Fenêtre d'ouverture/téléchargement")
                ),
                div(style = "flex: 0 0 auto;",
                    downloadButton(
                      outputId = paste0("download_ab1_", i),
                      label = "📂 Ouvrir",
                      style = "background-color: #28a745; color: white; border: none; padding: 6px 12px; border-radius: 4px; font-size: 12px;",
                      title = paste0("Ouvrir ", ab1_name)
                    )
                )
            )
          )
        } else {
          # LOCAL : Bouton d'ouverture directe (existant)
          button_list[[i]] <- div(
            style = "margin-bottom: 8px; padding: 8px; background: #f8f9fa; border: 1px solid #28a745; border-radius: 4px;",

            div(style = "display: flex; justify-content: space-between; align-items: center;",
                div(style = "flex: 1;",
                    tags$strong(style = "color: #28a745;", "✅ ", ab1_name)
                ),
                div(style = "flex: 0 0 auto;",
                    actionButton(
                      inputId = paste0("open_ab1_", i),
                      label = "📂 Ouvrir",
                      style = "background-color: #28a745; color: white; border: none; padding: 6px 12px; border-radius: 4px; font-size: 12px;",
                      title = paste0("Ouvrir ", ab1_name, " avec l'application par défaut")
                    )
                )
            )
          )
        }
      } else {
        # Fichier manquant
        button_list[[i]] <- div(
          style = "margin-bottom: 8px; padding: 8px; background: #fff5f5; border: 1px solid #dc3545; border-radius: 4px;",

          div(style = "display: flex; justify-content: space-between; align-items: center;",
              div(style = "flex: 1;",
                  tags$strong(style = "color: #dc3545;", "❌ ", ab1_name)
              ),
              div(style = "flex: 0 0 auto;",
                  tags$button(
                    "❌ Indisponible",
                    style = "background-color: #6c757d; color: white; border: none; padding: 6px 12px; border-radius: 4px; font-size: 12px; cursor: not-allowed;",
                    disabled = TRUE
                  )
              )
          )
        )
      }
    }

    return(div(
      # Message d'aide simple
      div(style = "background: #e8f5e8; padding: 8px; border-radius: 4px; margin-bottom: 10px; font-size: 12px;",
          "💡 ", tags$strong("Sur le serveur :"), " Cliquez pour ouvrir la fenêtre de téléchargement. ",
          "Choisissez 'Ouvrir avec' pour ouvrir directement le fichier."),
      button_list
    ))
  })



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
    selected_files <- get_all_selected_files()

    if (length(selected_files) == 0) {
      showNotification("⚠️ Aucun fichier sélectionné pour l'alignement", type = "warning", duration = 3)
      return(NULL)
    }

    # Organiser les fichiers par fragments
    fragments_data <- list()
    for (file_path in selected_files) {
      fragment_type <- extract_fragment_type(file_path)
      if (is.null(fragment_type)) {
        fragment_type <- "unknown"
      }

      fragments_data[[length(fragments_data) + 1]] <- list(
        path = file_path,
        type = fragment_type
      )
    }

    # Trier par type de fragment
    fragments_data <- fragments_data[order(sapply(fragments_data, function(x) {
      type_order <- c("5p", "int", "int1", "int2", "int3", "int4", "int5", "3p", "unknown")
      match(x$type, type_order)
    }))]

    # Charger les séquences
    sequences <- list()
    for (i in seq_along(fragments_data)) {
      fragment <- fragments_data[[i]]

      # Lire la séquence
      lines <- readLines(fragment$path, warn = FALSE)
      seq_raw <- paste(lines, collapse = "")
      seq_clean <- clean_sequence(seq_raw)

      # Appliquer le reverse complément pour les fragments 3p
      if (fragment$type == "3p") {
        seq_clean <- reverse_complement(seq_clean)
      }

      sequences[[i]] <- Biostrings::DNAString(seq_clean)
    }

    # Stocker les informations de fragments pour utilisation ultérieure
    attr(sequences, "fragments_info") <- fragments_data

    return(sequences)
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
    selected_files <- get_all_selected_files()

    if (length(selected_files) > 0 && !is.null(seqs())) {
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

      # Calcul des informations d'alignement pour la visualisation globale
      alignments_info <- calculate_alignments_info(seqs(), data_xdna$seq)

      # Génération de la visualisation globale
      overview_viz <- generate_alignment_overview(
        sequence_length = length(data_xdna$seq),
        features_lines = data_xdna$features,
        restriction_sites_list = restriction_sites(),
        alignments_info = alignments_info,
        selected_files = get_all_selected_files()
      )

      # Calcul de la région d'affichage si demandé
      display_region <- if (input$show_restriction_context) {
        calculate_restriction_display_region(restriction_sites(), length(data_xdna$seq))
      } else {
        list(start = 1, end = length(data_xdna$seq))
      }

      # Génération de la légende des couleurs
      legend_content <- generate_color_legend(data_xdna$features, restriction_sites())

      # Ajout d'information sur la région affichée
      if (input$show_restriction_context && length(restriction_sites()) > 0) {
        region_info <- paste0(
          "<div style='background: #e3f2fd; padding: 8px; border-radius: 4px; margin-bottom: 10px; font-size: 12px;'>",
          "📍 <strong>Région affichée:</strong> ", display_region$start, " - ", display_region$end,
          " (", display_region$end - display_region$start + 1, " nt)",
          "<br>💡 Centrée sur les sites de restriction avec ±200nt de contexte",
          "</div>"
        )
        legend_content <- paste0(region_info, legend_content)
      } else if (input$show_restriction_context && length(restriction_sites()) == 0) {
        region_info <- paste0(
          "<div style='background: #fff3cd; padding: 8px; border-radius: 4px; margin-bottom: 10px; font-size: 12px;'>",
          "ℹ️ <strong>Aucun site de restriction trouvé.</strong> Affichage de la séquence complète.",
          "</div>"
        )
        legend_content <- paste0(region_info, legend_content)
      } else {
        region_info <- paste0(
          "<div style='background: #f8f9fa; padding: 8px; border-radius: 4px; margin-bottom: 10px; font-size: 12px;'>",
          "📍 <strong>Mode:</strong> Affichage de la séquence complète de référence",
          "</div>"
        )
        legend_content <- paste0(region_info, legend_content)
      }

      align_output <- character()
      text_output <- character()

      selected_files_paths <- get_all_selected_files()

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

        # MODIFICATION : Utilisation de la région d'affichage
        if (input$show_restriction_context && length(restriction_sites()) > 0) {
          # Extraction de la sous-séquence dans la région d'intérêt
          region_seq <- data_xdna$seq[display_region$start:display_region$end]
          region_length <- length(region_seq)

          # Ajustement des positions d'alignement relatives à la région d'affichage
          relative_start <- max(1, aln_start - display_region$start + 1)
          relative_end <- min(region_length, aln_end - display_region$start + 1)

          # Création des séquences avec gaps pour la région
          full_pattern <- create_full_pattern_with_gaps(pat_aligned, relative_start, region_length)
          full_subject <- as.character(region_seq)
          full_annot <- create_full_annotation_with_spaces(annot_aligned, relative_start, region_length)

          # Cartes de couleurs et restrictions pour la région
          colors <- build_sequence_color_map(data_xdna$features, display_region$start, region_length)
          restriction_positions <- build_restriction_position_map(region_length, restriction_sites(), display_region$start)

          region_info_text <- paste0(" (région ", display_region$start, "-", display_region$end,
                                     ", alignement: ", aln_start, "-", aln_end, ")")

          # AJOUT CRUCIAL : Passer l'offset de début pour les positions absolues
          display_start_offset <- display_region$start
        } else {
          # Mode normal : séquence complète
          full_pattern <- create_full_pattern_with_gaps(pat_aligned, aln_start, length(data_xdna$seq))
          full_subject <- as.character(data_xdna$seq)
          full_annot <- create_full_annotation_with_spaces(annot_aligned, aln_start, length(data_xdna$seq))

          colors <- build_sequence_color_map(data_xdna$features, 1, length(data_xdna$seq))
          restriction_positions <- build_restriction_position_map(length(data_xdna$seq), restriction_sites())

          region_info_text <- paste0(" (région alignée: ", aln_start, "-", aln_end, ")")

          # Pas d'offset en mode normal
          display_start_offset <- 1
        }

        # Déterminer le type de fragment pour ce fichier
        file_fragment_type <- if (!is.null(selected_files_paths) && i <= length(selected_files_paths)) {
          extract_fragment_type(selected_files_paths[i])
        } else {
          NULL
        }

        # Si pas de type détecté, utiliser un fallback
        if (is.null(file_fragment_type)) {
          file_fragment_type <- "Seq"
        }

        # Nom d'affichage du fichier
        file_display_name <- if (!is.null(selected_files_paths) && i <= length(selected_files_paths)) {
          basename(selected_files_paths[i])
        } else {
          paste0("Fichier_", i)
        }

        # Génération de l'alignement coloré AVEC le type de fragment
        blast_alignment <- generate_colored_alignment(
          full_pattern, full_subject, full_annot,
          aln_start,  # IMPORTANT : Toujours utiliser la position absolue d'origine
          colors,
          paste0(file_display_name, region_info_text),
          restriction_positions,
          display_start_offset,  # NOUVEAU PARAMÈTRE pour calculer les positions absolues
          file_fragment_type     # NOUVEAU PARAMÈTRE pour le type de fragment
        )

        align_output <- c(align_output, blast_alignment$html)
        text_output <- c(text_output, blast_alignment$text)
      }

      incProgress(0.2, detail = "Finalisation...")

      # Sauvegarde pour les exports
      alignment_data$results <- align_output
      alignment_data$text_version <- paste(text_output, collapse = "")

      # Structure HTML finale avec visualisation globale + alignements détaillés
      tags$div(
        # Visualisation globale en haut
        HTML(overview_viz),

        # Alignements détaillés en bas
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

      selected_files <- get_all_selected_files()
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

      selected_files <- get_all_selected_files()
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
