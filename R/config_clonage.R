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
    # Chemins pour développement local (Windows)
    config <- list(
      xdna_dir = "R:/Production/Labo YEAST/Demandes du service/carte_nouveaux_clonages",  # NOUVEAU CHEMIN WINDOWS
      seq_dir = "P:/SEQ",                 # Dossier racine séquences (.seq)
      environment = "development"
    )

    # Chemins alternatifs pour dev Linux/Mac
    if (.Platform$OS.type != "windows") {
      config$xdna_dir <- "/mnt/carte_nouveaux_clonages"  # NOUVEAU CHEMIN LINUX
      if (dir.exists("data/genbank")) {
        config$seq_dir <- "data/seq"
      }
    }

  } else {
    # ========== ENVIRONNEMENT DE PRODUCTION ==========
    # Chemins pour production (Docker/ShinyProxy)
    config <- list(
      xdna_dir = "/mnt/carte_nouveaux_clonages",      # NOUVEAU CHEMIN MONTÉ
      seq_dir = "/data/production/SEQ",               # Dossiers avec fichiers .seq
      environment = "production"
    )

    # Vérification des chemins de production et fallbacks
    if (!dir.exists(config$xdna_dir)) {
      # Essayer des chemins alternatifs
      if (dir.exists("/data/SEQ/carte_nouveaux_clonages")) {
        config$xdna_dir <- "/data/SEQ/carte_nouveaux_clonages"
        config$seq_dir <- "/data/SEQ"
      } else if (dir.exists("../data/production/SEQ/carte_nouveaux_clonages")) {
        config$xdna_dir <- "../data/production/SEQ/carte_nouveaux_clonages"
        config$seq_dir <- "../data/production/SEQ"
      } else if (dir.exists("../data/SEQ/carte_nouveaux_clonages")) {
        config$xdna_dir <- "../data/SEQ/carte_nouveaux_clonages"
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

  # Debug détaillé des montages - MISE À JOUR
  cat("\n🔍 Debug montages carte_nouveaux_clonages:\n")
  cat("   - /mnt existe:", dir.exists("/mnt"), "\n")
  cat("   - /mnt/carte_nouveaux_clonages existe:", dir.exists("/mnt/carte_nouveaux_clonages"), "\n")
  cat("   - R:/Production/Labo YEAST/... existe:", dir.exists("R:/Production/Labo YEAST/Demandes du service/carte_nouveaux_clonages"), "\n")
  cat("   - /data existe:", dir.exists("/data"), "\n")
  cat("   - /data/production existe:", dir.exists("/data/production"), "\n")
  cat("   - /data/production/SEQ existe:", dir.exists("/data/production/SEQ"), "\n")

  # Lister le contenu de /data si accessible
  if (dir.exists("/data")) {
    data_content <- list.dirs("/data", recursive = FALSE, full.names = FALSE)
    cat("   - Contenu de /data:", paste(data_content, collapse = ", "), "\n")
  }

  # Lister le contenu de /mnt si accessible
  if (dir.exists("/mnt")) {
    mnt_content <- list.dirs("/mnt", recursive = FALSE, full.names = FALSE)
    cat("   - Contenu de /mnt:", paste(mnt_content, collapse = ", "), "\n")
  }

  # Lister les points de montage
  cat("\n💽 Points de montage disponibles:\n")
  mount_points <- c("/", "/data", "/srv", "/tmp", "/var", "/home", "/mnt")
  for (mp in mount_points) {
    if (dir.exists(mp)) {
      content <- try(list.dirs(mp, recursive = FALSE, full.names = FALSE), silent = TRUE)
      if (!inherits(content, "try-error") && length(content) > 0) {
        cat("   -", mp, ":", paste(head(content, 5), collapse = ", "), "\n")
      }
    }
  }

  # Affichage du contenu des dossiers GenBank (si accessibles)
  if (dir.exists(config$xdna_dir)) {
    gb_files <- list.files(config$xdna_dir, pattern = "\\.gb$")
    cat("\n📄 Fichiers GenBank dans", config$xdna_dir, ":\n")
    cat("   - Fichiers .gb trouvés:", length(gb_files), "\n")
    if (length(gb_files) > 0) {
      cat("     Exemples:", paste(head(gb_files, 3), collapse = ", "), "\n")
    }
  } else {
    cat("\n❌ Dossier GenBank non accessible:", config$xdna_dir, "\n")
  }

  # Vérifications spécifiques pour les nouveaux chemins
  cat("\n🔍 Vérifications spécifiques nouveaux chemins:\n")

  # Chemin Windows
  windows_path <- "R:/Production/Labo YEAST/Demandes du service/carte_nouveaux_clonages"
  if (dir.exists(windows_path)) {
    gb_files_win <- list.files(windows_path, pattern = "\\.gb$")
    cat("   - Windows path accessible:", length(gb_files_win), "fichiers .gb\n")
    if (length(gb_files_win) > 0) {
      cat("     Exemples:", paste(head(gb_files_win, 3), collapse = ", "), "\n")
    }
  } else {
    cat("   - Windows path non accessible:", windows_path, "\n")
  }

  # Chemin Linux monté
  linux_path <- "/mnt/carte_nouveaux_clonages"
  if (dir.exists(linux_path)) {
    gb_files_linux <- list.files(linux_path, pattern = "\\.gb$")
    cat("   - Linux mount accessible:", length(gb_files_linux), "fichiers .gb\n")
    if (length(gb_files_linux) > 0) {
      cat("     Exemples:", paste(head(gb_files_linux, 3), collapse = ", "), "\n")
    }
  } else {
    cat("   - Linux mount non accessible:", linux_path, "\n")
  }

  # Affichage du contenu des dossiers séquences (si accessibles)
  if (dir.exists(config$seq_dir)) {
    seq_folders <- list.dirs(config$seq_dir, recursive = FALSE, full.names = FALSE)
    cat("\n📁 Dossiers séquences dans", config$seq_dir, ":\n")
    cat("   - Dossiers trouvés:", length(seq_folders), "\n")
    if (length(seq_folders) > 0) {
      cat("     Exemples:", paste(head(seq_folders, 3), collapse = ", "), "\n")
    }
  } else {
    cat("\n❌ Dossier séquences non accessible:", config$seq_dir, "\n")
  }

  cat("══════════════════════════════\n")
}
