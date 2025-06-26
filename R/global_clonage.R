library(shiny)
library(Biostrings)
library(ape)
library(shinythemes)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "C:/Users/cantoine/AppData/Local/R/cache/R/reticulate/uv/cache/archive-v0/EPPRbPwPsqpfrrt7Mtpv1/Scripts/python.exe")

# Import Python module
snapgene <- import("snapgene_reader")

# Dossiers fichiers
xdna_dir <- "P:/SEQ/Atest_cae"
seq_dir <- "P:/SEQ/Atest_cae"

# Fonction nettoyage séquence : supprime espaces, retours à la ligne et met en majuscule
clean_seq <- function(x) {
  x <- gsub("[\r\n ]", "", x)
  toupper(x)
}

# Fonction pour annoter mutations (barre verticale à la différence)
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

# Fonction pour annoter vecteurs (ici simple ligne de points)
annotate_vectors <- function(seq_length, features = NULL) {
  vec_line <- rep(".", seq_length)
  # Ajoute traitement features si nécessaire
  paste(vec_line, collapse = "")
}



