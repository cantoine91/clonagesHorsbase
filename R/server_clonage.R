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
      # Alignement local Smith-Waterman
      aln <- Biostrings::pairwiseAlignment(
        pattern = seqs()[[i]],
        subject = data_xdna$seq,
        type = "local",
        substitutionMatrix = nucleotideSubstitutionMatrix(match = 2, mismatch = -1),
        gapOpening = -2,
        gapExtension = -1
      )

      # Position de l'alignement dans la séquence complète
      aln_start <- start(subject(aln))
      aln_end <- end(subject(aln))

      # Séquences alignées localement
      pat_aligned <- as.character(pattern(aln))
      sub_aligned <- as.character(subject(aln))
      annot_aligned <- annotate_mutations(pat_aligned, sub_aligned)

      # Créer la séquence pattern complète (gaps partout sauf dans la région alignée)
      full_pattern <- create_full_pattern_sequence(pat_aligned, aln_start, length(data_xdna$seq))

      # Séquence de référence complète
      full_subject <- as.character(data_xdna$seq)

      # Annotation complète (espaces partout sauf dans la région alignée)
      full_annot <- create_full_annotation_sequence(annot_aligned, aln_start, length(data_xdna$seq))

      # Couleurs pour toute la séquence de référence (positions 1 à longueur totale)
      colors <- get_color_map(data_xdna$features, 1, length(data_xdna$seq))

      # Générer l'alignement style SerialCloner
      blast_alignment <- generate_blast_style_alignment(
        full_pattern, full_subject, full_annot, aln_start, colors,
        paste0(input$seq_files[i], " (région alignée: ", aln_start, "-", aln_end, ")")
      )

      align_output <- c(align_output, blast_alignment$html)
      text_output <- c(text_output, blast_alignment$text)
    }

    # Stocker les résultats pour les exports
    alignment_data$results <- align_output
    alignment_data$text_version <- paste(text_output, collapse = "")

    # Structure HTML avec légende fixe + conteneur scrollable
    tags$div(
      id = "align_results",
      HTML(legend_html),
      tags$div(
        class = "alignments-container",
        HTML(paste(align_output, collapse = ""))
      )
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
