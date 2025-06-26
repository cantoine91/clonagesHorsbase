library(shiny)
library(Biostrings)
library(reticulate)

xdna_dir <- "P:/SEQ/Atest_cae"
seq_dir <- "P:/SEQ/Atest_cae"

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

# Parser features GenBank basique, extrait positions et couleur (exemple fixe)
parse_features <- function(features_lines) {
  feats <- list()
  for (line in features_lines) {
    pos_match <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
    if (length(pos_match) > 0) {
      feats[[length(feats) + 1]] <- list(pos_raw = pos_match, color = sample(c("red","blue","green","orange","purple"), 1))
    }
  }
  feats
}

# Construire vecteur couleurs pour la séquence
get_color_map <- function(features_lines, aln_start, seq_length) {
  feats <- parse_features(features_lines)
  color_map <- rep(NA, seq_length)
  for (feat in feats) {
    pos_str <- feat$pos_raw
    if (is.null(pos_str)) next
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

# Créer ligne HTML caractère par caractère monospace, couleur + tooltip
htmlify_line <- function(line, colors = NULL, start_pos = NULL) {
  chars <- strsplit(line, "")[[1]]
  spans <- character(length(chars))
  for (i in seq_along(chars)) {
    style <- "font-family: Courier New, monospace;"
    if (!is.null(colors) && !is.na(colors[i])) {
      style <- paste0(style, "color:", colors[i], ";")
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

      aln_start <- start(subject(aln))

      colors <- get_color_map(data_xdna$features, aln_start, nchar(sub))

      pattern_html <- htmlify_line(pat)
      annot_html <- htmlify_line(annot)
      subject_colored <- htmlify_line(sub, colors, aln_start)

      align_block <- paste0(
        "<b>Alignement avec ", input$seq_files[i], " :</b>\n",
        pattern_html, "\n",
        annot_html, "\n",
        subject_colored, "\n",
        "<i>Score : ", score(aln), "</i>"
      )

      align_output <- c(align_output, align_block)
    }

    tags$div(
      id = "align_results",
      style = "font-family: 'Courier New', monospace; background: #f8f8f8; border: 1px solid #ddd; padding: 10px; white-space: pre; overflow-x: auto; max-width: 100%; max-height: 400px;",
      HTML(paste(align_output, collapse = "\n\n"))
    )
  })

}
