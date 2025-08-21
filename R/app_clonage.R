# ==============================================================================
# APP_CLONAGE.R - Version nettoyée
# Point d'entrée principal de l'application HGX
# Charge tous les modules et lance l'application Shiny
# ==============================================================================

#' Application HGX - Analyseur de séquences de clonage
#'
#' Cette application permet d'analyser et comparer des séquences de clonage
#' avec des cartes de référence GenBank. Elle offre les fonctionnalités suivantes :
#'
#' FONCTIONNALITÉS PRINCIPALES :
#' - Chargement de cartes GenBank de référence (.gb/.txt)
#' - Recherche automatique de fichiers de séquences (.seq)
#' - Alignements locaux et globaux avec visualisation colorée
#' - Détection et affichage des sites de restriction
#' - Gestion des fragments 5', internes et 3'
#' - Export des résultats en HTML et FASTA
#' - Ouverture automatique des chromatogrammes (.ab1)
#'
#' ARCHITECTURE :
#' - config_clonage.R  : Configuration automatique des environnements
#' - global_clonage.R  : Fonctions utilitaires et algorithmes
#' - ui_clonage.R      : Interface utilisateur
#' - server_clonage.R  : Logique serveur et traitement des données
#'
#' UTILISATION :
#' 1. Sélectionner une carte GenBank de référence
#' 2. Rechercher des fichiers .seq par nom de plaque et mot-clé
#' 3. Choisir les enzymes de restriction à analyser
#' 4. Lancer l'alignement et visualiser les résultats
#' 5. Exporter les rapports selon les besoins

# Chargement des modules dans l'ordre de dépendance
source("config_clonage.R")    # Configuration des chemins
source("global_clonage.R")    # Fonctions utilitaires
source("ui_clonage.R")        # Interface utilisateur
source("server_clonage.R")    # Logique serveur

# Lancement de l'application Shiny
shinyApp(ui = ui_clonage, server = server_clonage)
