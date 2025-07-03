# ==============================================================================
# CONFIG_CLONAGE.R
# Configuration automatique des environnements pour HGX Clonage
# Détection automatique dev/production et chemins appropriés
# ==============================================================================

#' Détection automatique de l'environnement d'exécution
#' @return Chaîne indiquant l'environnement ("development" ou "production")
detect_environment <- function() {
  # Méthode 1: Vérifier si on est dans un conteneur Docker
  if (file.exists("/.dockerenv")) {
    return("production")
  }

  # Méthode 2: Vérifier les chemins Windows (dev local)
  if (.Platform$OS.type == "windows") {
    return("development")
  }

  # Méthode 3: Vérifier si le répertoire de production existe
  if (dir.exists("/data/production/SEQ")) {
    return("production")
  }

  # Méthode 4: Vérifier si on est dans un environnement ShinyProxy/serveur
  if (dir.exists("/srv/shiny-server")) {
    return("production")
  }

  # Par défaut: développement
  return("development")
}

#' Configuration des chemins selon l'environnement détecté
#' @return Liste avec chemins configurés et informations d'environnement
get_config <- function() {
  env <- detect_environment()

  if (env == "development") {
    # ========== ENVIRONNEMENT DE DÉVELOPPEMENT ==========
    # Chemins pour développement local (Windows/Mac/Linux)
    config <- list(
      xdna_dir = "P:/SEQ/Atest_cae",      # Dossier fichiers GenBank (.gb)
      seq_dir = "P:/SEQ",                 # Dossier racine séquences (.seq)
      environment = "development"
    )

    # Chemins alternatifs pour dev Linux/Mac
    if (.Platform$OS.type != "windows") {
      if (dir.exists("data/genbank")) {
        config$xdna_dir <- "data/genbank"
        config$seq_dir <- "data/seq"
      }
    }

  } else {
    # ========== ENVIRONNEMENT DE PRODUCTION ==========
    # Chemins pour production (Docker/ShinyProxy)
    config <- list(
      xdna_dir = "/data/production/SEQ/Atest_cae",    # Fichiers GenBank (.gb)
      seq_dir = "/data/production/SEQ",               # Dossiers avec fichiers .seq
      environment = "production"
    )

    # Vérification des chemins de production et fallbacks
    if (!dir.exists(config$xdna_dir)) {
      # Essayer des chemins alternatifs
      if (dir.exists("/data/SEQ/Atest_cae")) {
        config$xdna_dir <- "/data/SEQ/Atest_cae"
        config$seq_dir <- "/data/SEQ"
      } else if (dir.exists("../data/production/SEQ/Atest_cae")) {
        config$xdna_dir <- "../data/production/SEQ/Atest_cae"
        config$seq_dir <- "../data/production/SEQ"
      } else if (dir.exists("../data/SEQ/Atest_cae")) {
        config$xdna_dir <- "../data/SEQ/Atest_cae"
        config$seq_dir <- "../data/SEQ"
      }
    }
  }

  return(config)
}

#' Validation de la configuration
#' @param config Configuration à valider
#' @return TRUE si la configuration est valide, FALSE sinon
validate_config <- function(config) {
  xdna_exists <- dir.exists(config$xdna_dir)
  seq_exists <- dir.exists(config$seq_dir)

  if (!xdna_exists) {
    warning("⚠️ Dossier GenBank non trouvé: ", config$xdna_dir)
  }

  if (!seq_exists) {
    warning("⚠️ Dossier séquences non trouvé: ", config$seq_dir)
  }

  return(xdna_exists && seq_exists)
}

#' Affichage des informations de configuration (pour debug)
#' @param config Configuration à afficher
display_config_info <- function(config) {
  cat("🔧 Configuration HGX Clonage\n")
  cat("══════════════════════════════\n")
  cat("Environnement détecté:", config$environment, "\n")
  cat("Système d'exploitation:", .Platform$OS.type, "\n")
  cat("Docker detecté:", file.exists("/.dockerenv"), "\n")
  cat("\n📁 Chemins configurés:\n")
  cat("   - GenBank (.gb):", config$xdna_dir, "\n")
  cat("   - Séquences (.seq):", config$seq_dir, "\n")
  cat("\n✅ Validation:\n")
  cat("   - GenBank accessible:", dir.exists(config$xdna_dir), "\n")
  cat("   - Séquences accessibles:", dir.exists(config$seq_dir), "\n")

  # Affichage du contenu des dossiers (si accessibles)
  if (dir.exists(config$xdna_dir)) {
    gb_files <- list.files(config$xdna_dir, pattern = "\\.gb$")
    cat("   - Fichiers .gb trouvés:", length(gb_files), "\n")
    if (length(gb_files) > 0) {
      cat("     Exemples:", paste(head(gb_files, 3), collapse = ", "), "\n")
    }
  }

  if (dir.exists(config$seq_dir)) {
    seq_folders <- list.dirs(config$seq_dir, recursive = FALSE, full.names = FALSE)
    cat("   - Dossiers séquences:", length(seq_folders), "\n")
    if (length(seq_folders) > 0) {
      cat("     Exemples:", paste(head(seq_folders, 3), collapse = ", "), "\n")
    }
  }
  cat("══════════════════════════════\n")
}
