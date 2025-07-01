# config.R - Configuration automatique des environnements
# À ajouter au début de global_clonage.R ou dans un fichier séparé

# Fonction pour détecter l'environnement
detect_environment <- function() {
  # Méthode 1: Vérifier si on est dans un conteneur Docker
  if (file.exists("/.dockerenv")) {
    return("production")
  }

  # Méthode 2: Vérifier les chemins Windows (dev local)
  if (.Platform$OS.type == "windows") {
    return("development")
  }

  # Méthode 3: Vérifier si le répertoire de prod existe
  if (dir.exists("/srv/shiny-server")) {
    return("production")
  }

  # Par défaut: développement
  return("development")
}

# Configuration des chemins selon l'environnement
get_config <- function() {
  env <- detect_environment()

  if (env == "development") {
    # Chemins pour développement local (Windows/Mac/Linux)
    list(
      xdna_dir = "P:/SEQ/Atest_cae",  # ou "data/genbank" si vous avez des données locales
      seq_dir = "P:/SEQ/Atest_cae",   # ou "data/seq"
      environment = "development"
    )
  } else {
    # Chemins pour production (Docker)
    list(
      xdna_dir = "../data/genbank",
      seq_dir = "../data/seq",
      environment = "production"
    )
  }
}

# Utilisation
config <- get_config()
xdna_dir <- config$xdna_dir
seq_dir <- config$seq_dir

# Debug: afficher l'environnement détecté
cat("Environnement détecté:", config$environment, "\n")
cat("xdna_dir:", xdna_dir, "\n")
cat("seq_dir:", seq_dir, "\n")
