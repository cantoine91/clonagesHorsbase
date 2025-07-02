library(shiny)
library(shinythemes)

ui_clonage <- shinyUI(navbarPage(
  title = div(style = "color: white; font-weight: bold; font-size: 20px;", "HGX"),
  id = "navbar",
  theme = shinytheme("united"),
  windowTitle = "Titre fenêtre navigateur",

  # Ajouter les styles pour l'affichage côte à côte et les tooltips
  header = tags$head(
    tags$style(HTML("
      .results-container {
        display: flex;
        gap: 20px;
        align-items: flex-start;
      }

      .legend-container {
        flex: 0 0 300px;
        background: #f8f9fa;
        padding: 15px;
        border-radius: 5px;
        border: 1px solid #dee2e6;
        font-family: Arial, sans-serif;
        font-size: 13px;
        max-height: 600px;
        overflow-y: auto;
        position: sticky;
        top: 20px;
      }

      .alignments-container {
        flex: 1;
        background: #f8f8f8;
        padding: 15px;
        border-radius: 5px;
        border: 1px solid #ddd;
        overflow-x: auto;
        max-height: 600px;
        font-family: 'Courier New', monospace;
        line-height: 1.4;
      }

      .alignment-block {
        margin-bottom: 25px;
        padding: 15px;
        background: #ffffff;
        border: 1px solid #e0e0e0;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        white-space: pre;
        overflow-x: auto;
        min-width: fit-content;
      }

      .alignment-title {
        color: #2c3e50;
        font-weight: bold;
        margin-bottom: 10px;
        font-family: Arial, sans-serif;
        border-bottom: 2px solid #b22222;
        padding-bottom: 5px;
        white-space: normal;
      }

      /* Styles pour les tooltips */
      .nucleotide-cell {
        position: relative;
        cursor: pointer;
      }

      .nucleotide-cell:hover {
        background-color: #e3f2fd !important;
        border: 1px solid #2196f3 !important;
        z-index: 1000;
        border-radius: 3px !important;
      }

      .tooltip {
        visibility: hidden;
        background-color: #333;
        color: white;
        text-align: center;
        border-radius: 4px;
        padding: 5px 8px;
        position: absolute;
        z-index: 1001;
        bottom: 125%;
        left: 50%;
        margin-left: -60px;
        opacity: 0;
        transition: opacity 0.3s;
        font-size: 11px;
        white-space: nowrap;
        pointer-events: none;
        box-shadow: 0 2px 8px rgba(0,0,0,0.3);
        min-width: 120px;
      }

      .tooltip::after {
        content: '';
        position: absolute;
        top: 100%;
        left: 50%;
        margin-left: -5px;
        border-width: 5px;
        border-style: solid;
        border-color: #333 transparent transparent transparent;
      }

      .nucleotide-cell:hover .tooltip {
        visibility: visible;
        opacity: 1;
      }
    ")),

    # JavaScript pour la fonctionnalité de copie
    tags$script(HTML("
      function copyToClipboard() {
        const alignResults = document.querySelector('.alignments-container');
        if (!alignResults) return;

        const textContent = alignResults.innerText || alignResults.textContent;

        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(textContent).then(function() {
            alert('Alignement copié dans le presse-papiers !');
          }).catch(function(err) {
            console.error('Erreur lors de la copie: ', err);
          });
        }
      }
    "))
  ),

  tabPanel("Clonages Hors Base",

           # Interface simple et entièrement adaptative
           div(style = "width: 100%; margin: 20px auto; padding: 20px;",

               # Fichiers d'entrée
               wellPanel(
                 h4("📁 Fichiers d'entrée", style = "color: #b22222; margin-top: 0;"),
                 fluidRow(
                   column(6,
                          actionButton("refresh_files", "🔄", style = "float: right; margin-bottom: 5px; background: transparent; border: 1px solid #ccc;"),
                          selectInput("carte_xdna", "Carte de référence (.gb):", choices = NULL)
                   ),
                   column(6,
                          selectInput("seq_files", "Séquences test (.seq):", choices = NULL, multiple = TRUE)
                   )
                 )
               ),

               # Enzymes de restriction
               wellPanel(
                 h4("🧬 Enzymes de restriction", style = "color: #b22222; margin-top: 0;"),
                 fluidRow(
                   column(6,
                          selectInput("enzyme1", "Enzyme 1:",
                                      choices = c("Aucune" = "", names(get_restriction_enzymes())),
                                      selected = "")
                   ),
                   column(6,
                          selectInput("enzyme2", "Enzyme 2:",
                                      choices = c("Aucune" = "", names(get_restriction_enzymes())),
                                      selected = "")
                   )
                 ),
                 div(style = "background: #e8f4f8; padding: 8px; border-radius: 4px; margin-top: 10px; font-family: monospace; font-size: 12px;",
                     textOutput("restriction_info")
                 )
               ),

               # Bouton d'alignement
               div(style = "text-align: center; margin: 20px 0;",
                   actionButton("align_btn", "🔬 Lancer l'alignement",
                                style = "background-color: #b22222; color: white; font-size: 16px; padding: 12px 30px; border: none; border-radius: 5px; font-weight: bold;")
               ),

               # Informations
               wellPanel(
                 h4("📊 Informations", style = "color: #b22222; margin-top: 0;"),
                 verbatimTextOutput("seq_xdna"),
                 verbatimTextOutput("seqs_selected")
               ),

               # Instructions pour les tooltips
               conditionalPanel(
                 condition = "output.align_results",
                 div(style = "background: #e8f5e8; padding: 10px; border-radius: 4px; margin: 10px 0; border-left: 4px solid #4caf50;",
                     tags$p(style = "margin: 0; font-size: 13px; color: #2e7d32;",
                            "💡 ", tags$strong("Astuce:"), " Survolez les nucléotides dans l'alignement pour voir leur position exacte dans la séquence !")
                 )
               ),

               # Boutons d'export
               conditionalPanel(
                 condition = "output.align_results",
                 div(style = "margin: 20px 0; text-align: center;",
                     tags$button("📋 Copier", onclick = "copyToClipboard()",
                                 style = "background: #27ae60; color: white; border: none; padding: 8px 15px; margin-right: 10px; border-radius: 4px;"),
                     downloadButton("download_txt", "💾 Télécharger TXT",
                                    style = "background: #2c3e50; color: white; border: none; padding: 8px 15px; border-radius: 4px;")
                 )
               ),

               # Résultats d'alignement
               uiOutput("align_results")
           )
  ),

  tabPanel("TEST",
           h4("Contenu ici...")
  )
))
