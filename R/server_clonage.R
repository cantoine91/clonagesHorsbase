library(shiny)
library(Biostrings)
library(reticulate)

xdna_dir <- "P:/SEQ/Atest_cae"
seq_dir <- "P:/SEQ/Atest_cae"

server_clonage <- function(input, output, session) {
  data_xdna <- reactiveValues(seq = NULL, features = NULL)
  alignment_data <- reactiveValues(results = NULL, text_version = NULL)

  observeEvent(input$align_btn, {
    req(input$carte_xdna)
    fichier <- file.path(xdna_dir, input$carte_xdna)
    gb_lines <- readLines(fichier, warn = FALSE)

    origin_line <- grep("^ORIGIN", gb_lines, ignore.case = TRUE)

    # --- Features
    features_block <- gb_lines[1:(origin_line - 1)]
    features_lines <- features_block[grep("^\\s{5}|^\\s{21}", features_block)]
    data_xdna$features <- features_lines

    # --- Séquence propre
    seq_lines <- gb_lines[(origin_line + 1):length(gb_lines)]
    seq_raw <- paste(seq_lines, collapse = "")
    seq_clean <- gsub("[^acgtACGTnN]", "", seq_raw)
    data_xdna$seq <- Biostrings::DNAString(toupper(seq_clean))
  })

  seqs <- eventReactive(input$align_btn, {
    req(input$seq_files)
    paths <- file.path(seq_dir, input$seq_files)
    lapply(paths, function(f) {
      lines <- readLines(f, warn = FALSE)
      seq_raw <- paste(lines, collapse = "")
      seq_clean <- clean_seq(seq_raw)
      Biostrings::DNAString(seq_clean)
    })
  })

  output$seq_xdna <- renderPrint({
    req(data_xdna$seq)
    cat("Séquence GenBank :\n")
    print(data_xdna$seq)
    cat("\nAnnotations (Features) :\n")
    print(data_xdna$features)
  })

  output$seqs_selected <- renderPrint({
    req(seqs())
    cat("Séquences chargées des fichiers .seq :\n")
    for (i in seq_along(seqs())) {
      cat(" - ", input$seq_files[i], ": ", as.character(seqs()[[i]]), "\n")
    }
  })

  output$align_results <- renderUI({
    req(data_xdna$seq, seqs())

    # Générer la légende des couleurs
    legend_html <- generate_legend(data_xdna$features)

    align_output <- character()
    text_output <- character()

    for (i in seq_along(seqs())) {
      aln <- Biostrings::pairwiseAlignment(pattern = seqs()[[i]], subject = data_xdna$seq)

      pat <- as.character(pattern(aln))
      sub <- as.character(subject(aln))
      annot <- annotate_mutations(pat, sub)

      aln_start <- start(subject(aln))

      colors <- get_color_map(data_xdna$features, aln_start, nchar(sub))

      # Générer la règle de numérotation
      ruler <- generate_ruler(aln_start, nchar(sub))

      pattern_html <- htmlify_line(pat)
      annot_html <- htmlify_line(annot)
      subject_colored <- htmlify_line(sub, colors, aln_start)
      ruler_html <- paste0("<span class='ruler-line'>", htmlify_line(ruler), "</span>")

      align_block <- paste0(
        "<div class='alignment-block'>",
        "<div class='alignment-title'>Alignement avec ", input$seq_files[i], "</div>",
        ruler_html, "\n",
        pattern_html, "\n",
        annot_html, "\n",
        subject_colored, "\n",
        "<div class='alignment-score'>Score : ", score(aln), "</div>",
        "</div>"
      )

      # Version texte pour export
      text_block <- paste0(
        "=== Alignement avec ", input$seq_files[i], " ===\n",
        "Score: ", score(aln), "\n",
        "Position: ", aln_start, "-", aln_start + nchar(sub) - 1, "\n\n",
        "Ruler:   ", ruler, "\n",
        "Pattern: ", pat, "\n",
        "Match:   ", annot, "\n",
        "Subject: ", sub, "\n\n"
      )

      align_output <- c(align_output, align_block)
      text_output <- c(text_output, text_block)
    }

    # Stocker les résultats pour les exports
    alignment_data$results <- align_output
    alignment_data$text_version <- paste(text_output, collapse = "")

    tags$div(
      id = "align_results",
      HTML(paste0(legend_html, paste(align_output, collapse = "")))
    )
  })

  # Download handlers
  output$download_txt <- downloadHandler(
    filename = function() {
      paste0("alignement_", Sys.Date(), ".txt")
    },
    content = function(file) {
      req(alignment_data$text_version)

      # En-tête
      header <- paste0(
        "=== RESULTATS D'ALIGNEMENT - HGX ===\n",
        "Date: ", Sys.time(), "\n",
        "Carte GenBank: ", input$carte_xdna, "\n",
        "Fichiers séquences: ", paste(input$seq_files, collapse = ", "), "\n",
        "Longueur séquence référence: ", length(data_xdna$seq), " nt\n\n",
        "LEGENDE DES COULEURS:\n"
      )

      # Légende en format texte
      feats <- parse_features(data_xdna$features)
      legend_text <- ""
      for (feat in feats) {
        if (!is.na(feat$pos_raw)) {
          bounds <- as.numeric(unlist(strsplit(feat$pos_raw, "\\.\\.")))
          if (length(bounds) == 2) {
            name_display <- if (feat$name != "") feat$name else "Feature sans nom"
            legend_text <- paste0(legend_text, "- ", name_display, " (", bounds[1], "-", bounds[2], ") - Couleur: ", feat$color, "\n")
          }
        }
      }

      content <- paste0(header, legend_text, "\n", alignment_data$text_version)
      writeLines(content, file, useBytes = TRUE)
    }
  )

  output$download_html <- downloadHandler(
    filename = function() {
      paste0("alignement_", Sys.Date(), ".html")
    },
    content = function(file) {
      req(alignment_data$results)

      legend_html <- generate_legend(data_xdna$features)

      html_content <- paste0(
        "<!DOCTYPE html>\n<html>\n<head>\n",
        "<title>Résultats d'alignement - HGX</title>\n",
        "<meta charset='UTF-8'>\n",
        "<style>\n",
        "body { font-family: Arial, sans-serif; margin: 20px; }\n",
        "h1 { color: #b22222; }\n",
        ".info { background: #f0f0f0; padding: 10px; margin-bottom: 15px; border-radius: 5px; }\n",
        ".alignment-results { font-family: 'Courier New', monospace; background: #f8f8f8; padding: 15px; border-radius: 5px; }\n",
        ".alignment-block { margin-bottom: 25px; padding: 15px; background: white; border: 1px solid #ddd; border-radius: 5px; }\n",
        ".alignment-title { color: #2c3e50; font-weight: bold; margin-bottom: 10px; border-bottom: 2px solid #3498db; padding-bottom: 5px; }\n",
        ".ruler-line { color: #666; background: #f9f9f9; }\n",
        "</style>\n</head>\n<body>\n",
        "<h1>Résultats d'alignement - HGX</h1>\n",
        "<div class='info'>\n",
        "<strong>Date:</strong> ", Sys.time(), "<br>\n",
        "<strong>Carte GenBank:</strong> ", input$carte_xdna, "<br>\n",
        "<strong>Fichiers séquences:</strong> ", paste(input$seq_files, collapse = ", "), "<br>\n",
        "<strong>Longueur séquence référence:</strong> ", length(data_xdna$seq), " nt\n",
        "</div>\n",
        "<div class='alignment-results'>\n",
        legend_html,
        paste(alignment_data$results, collapse = ""),
        "</div>\n</body>\n</html>"
      )

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

      # Séquence référence
      fasta_content <- c(fasta_content,
                         paste0(">Sequence_Reference_", gsub("\\.gb$", "", input$carte_xdna)),
                         as.character(data_xdna$seq))

      # Séquences alignées
      for (i in seq_along(seqs())) {
        seq_name <- gsub("\\.seq$", "", input$seq_files[i])
        fasta_content <- c(fasta_content,
                           paste0(">Sequence_", seq_name),
                           as.character(seqs()[[i]]))
      }

      writeLines(fasta_content, file)
    }
  )
}
