library(shiny)
library(Biostrings)
library(reticulate)

server_clonage <- function(input, output, session) {
  data_xdna <- reactiveValues(seq = NULL, features = NULL)

  observeEvent(input$align_btn, {
    req(input$carte_xdna)
    fichier <- file.path(xdna_dir, input$carte_xdna)
    gb_lines <- readLines(fichier, warn = FALSE)
    origin_line <- grep("^ORIGIN", gb_lines, ignore.case = TRUE)
    if (length(origin_line) == 0) {
      showNotification("Mot-clé 'ORIGIN' non trouvé dans le fichier .gb", type = "error")
      return(NULL)
    }
    sequence_lines <- gb_lines[(origin_line + 1):length(gb_lines)]
    raw_seq <- paste(sequence_lines, collapse = "")
    sequence_text <- gsub("[^ACGTNacgtn]", "", raw_seq)
    sequence_text <- clean_seq(sequence_text)
    data_xdna$seq <- Biostrings::DNAString(sequence_text)
    data_xdna$features <- grep("^     ", gb_lines, value = TRUE)
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
    align_output <- character()

    for (i in seq_along(seqs())) {
      aln <- Biostrings::pairwiseAlignment(pattern = seqs()[[i]], subject = data_xdna$seq)

      pat <- as.character(pattern(aln))
      sub <- as.character(subject(aln))
      annot <- annotate_mutations(pat, sub)
      vec_line <- annotate_vectors(nchar(pat), data_xdna$features)

      align_block <- paste0(
        "Alignement avec ", input$seq_files[i], ":\n",
        pat, "\n",
        annot, "\n",
        sub, "\n",
        vec_line, "\n",
        "Score : ", score(aln)
      )

      align_output <- c(align_output, align_block)
    }

    HTML(paste0("<pre style='font-family: Courier New, monospace;'>",
                paste(align_output, collapse = "\n\n"),
                "</pre>"))
  })
}
