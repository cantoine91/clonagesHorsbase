library(shiny)
library(shinythemes)

xdna_dir <- "P:/SEQ/Atest_cae"
seq_dir <- "P:/SEQ/Atest_cae"

ui_clonage <- shinyUI(navbarPage(
  title = div(style = "color: white; font-weight: bold; font-size: 20px;", "HGX"),
  id = "navbar",
  theme = shinytheme("united"),
  windowTitle = "Titre fenêtre navigateur",
  header = tags$head(
    tags$style(HTML("
      .navbar { background-color: #b22222 !important; }
      .navbar .navbar-brand { color: white !important; }
      .navbar .navbar-nav > li > a { color: white !important; font-weight: bold; }
      input[type='text']:focus, textarea:focus {
        border: 2px solid #b22222 !important;
        box-shadow: 0 0 5px rgba(178,34,34,0.5);
        outline: none;
      }
      input[type='checkbox']:checked { accent-color: #b22222; }

      /* Full width container for inputs */
      #input_container {
        width: 100%;
        padding: 10px 15px;
        background: #f5f5f5;
        border-bottom: 1px solid #ddd;
        box-sizing: border-box;
        display: flex;
        flex-wrap: wrap;
        gap: 15px;
        align-items: center;
      }

      /* Inputs take equal width on bigger screens */
      #input_container > * {
        flex-grow: 1;
        min-width: 200px;
      }

      /* Export buttons container */
      #export_buttons {
        display: flex;
        gap: 10px;
        margin-bottom: 15px;
        padding: 10px;
        background: #f9f9f9;
        border-radius: 5px;
        border: 1px solid #ddd;
        flex-wrap: wrap;
      }

      .export-btn {
        background-color: #2c3e50 !important;
        color: white !important;
        border: none !important;
        padding: 8px 12px !important;
        border-radius: 4px !important;
        font-size: 12px !important;
        cursor: pointer !important;
        transition: background-color 0.3s !important;
      }

      .export-btn:hover {
        background-color: #34495e !important;
      }

      .copy-btn {
        background-color: #27ae60 !important;
      }

      .copy-btn:hover {
        background-color: #2ecc71 !important;
      }

      /* Full width results container */
      #results_container {
        width: 100%;
        padding: 15px;
        box-sizing: border-box;
      }

      /* Conteneur principal pour la légende + alignements */
      #align_results {
        width: 100% !important;
        max-width: 100% !important;
        border: 1px solid #ddd !important;
        border-radius: 5px !important;
        background: #f8f8f8 !important;
        position: relative !important;
        overflow: hidden !important;
      }

      /* Style pour la légende des couleurs - FIXE en haut */
      #legend-fixed {
        position: sticky !important;
        top: 0 !important;
        z-index: 1000 !important;
        margin: 0 !important;
        padding: 15px !important;
        background: #f8f9fa !important;
        border-bottom: 2px solid #dee2e6 !important;
        border-radius: 5px 5px 0 0 !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1) !important;
        font-family: Arial, sans-serif !important;
      }

      #legend-fixed b {
        color: #495057 !important;
        font-size: 14px !important;
      }

      /* Conteneur des alignements avec scroll horizontal */
      .alignments-container {
        padding: 15px !important;
        overflow-x: auto !important;
        overflow-y: auto !important;
        max-height: 600px !important;
        font-family: 'Courier New', monospace !important;
        line-height: 1.4 !important;
        width: 100% !important;
        box-sizing: border-box !important;
      }

      /* Style pour les alignements individuels */
      .alignment-block {
        margin-bottom: 25px !important;
        padding: 15px !important;
        background: #ffffff !important;
        border: 1px solid #e0e0e0 !important;
        border-radius: 5px !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1) !important;
        white-space: pre !important;
        overflow-x: auto !important;
        min-width: fit-content !important;
      }

      /* Assurer que chaque ligne d'alignement ne wrap pas */
      .alignment-block > div {
        white-space: nowrap !important;
        overflow-x: visible !important;
      }

      .alignment-title {
        color: #2c3e50 !important;
        font-weight: bold !important;
        margin-bottom: 10px !important;
        font-family: Arial, sans-serif !important;
        border-bottom: 2px solid #b22222 !important;
        padding-bottom: 5px !important;
        white-space: normal !important;
      }

      .alignment-score {
        color: #7f8c8d !important;
        font-style: italic !important;
        margin-top: 8px !important;
        font-family: Arial, sans-serif !important;
        white-space: normal !important;
      }

      /* Style pour la règle de numérotation */
      .ruler-line {
        color: #666 !important;
        font-weight: normal !important;
        border-bottom: 1px solid #ddd !important;
        margin-bottom: 2px !important;
        padding-bottom: 2px !important;
        background: transparent !important;
        display: block !important;
        font-family: 'Courier New', monospace !important;
        line-height: 1.4 !important;
      }

      /* Lignes de séquences alignées */
      .sequence-line {
        font-family: 'Courier New', monospace !important;
        line-height: 1.4 !important;
        margin: 0 !important;
        padding: 0 !important;
        white-space: nowrap !important;
        display: block !important;
      }

      /* Amélioration des tooltips */
      span[title]:hover {
        background-color: #ffffcc !important;
        border-radius: 2px !important;
      }

      /* Style pour les verbatimTextOutput */
      .shiny-text-output {
        background: #f9f9f9;
        border: 1px solid #ddd;
        border-radius: 5px;
        padding: 10px;
        margin-bottom: 15px;
        font-family: 'Courier New', monospace;
        font-size: 12px;
        max-height: 200px;
        overflow-y: auto;
      }

      /* Assurer que la police monospace est bien appliquée dans les alignements */
      .alignments-container * {
        font-family: 'Courier New', monospace !important;
      }

      .alignment-title, .alignment-score {
        font-family: Arial, sans-serif !important;
      }

      /* Forcer le contenu à ne pas déborder */
      .alignments-container span {
        display: inline !important;
        white-space: nowrap !important;
      }

      /* Style pour les notifications */
      .notification {
        position: fixed;
        top: 20px;
        right: 20px;
        padding: 10px 15px;
        background: #27ae60;
        color: white;
        border-radius: 5px;
        z-index: 9999;
        display: none;
      }

      /* Responsive design */
      @media (max-width: 768px) {
        .alignments-container {
          font-size: 10px !important;
        }
      }
    ")),

    # JavaScript pour la fonctionnalité de copie
    tags$script(HTML("
      function copyToClipboard() {
        const alignResults = document.querySelector('.alignments-container');
        if (!alignResults) return;

        // Créer une version texte propre
        const textContent = alignResults.innerText || alignResults.textContent;

        // Utiliser l'API moderne du clipboard
        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(textContent).then(function() {
            showNotification('Alignement copié dans le presse-papiers !');
          }).catch(function(err) {
            console.error('Erreur lors de la copie: ', err);
            fallbackCopy(textContent);
          });
        } else {
          fallbackCopy(textContent);
        }
      }

      function fallbackCopy(text) {
        const textArea = document.createElement('textarea');
        textArea.value = text;
        document.body.appendChild(textArea);
        textArea.select();
        try {
          document.execCommand('copy');
          showNotification('Alignement copié dans le presse-papiers !');
        } catch (err) {
          console.error('Erreur lors de la copie: ', err);
          showNotification('Erreur lors de la copie', 'error');
        }
        document.body.removeChild(textArea);
      }

      function showNotification(message, type = 'success') {
        let notification = document.querySelector('.notification');
        if (!notification) {
          notification = document.createElement('div');
          notification.className = 'notification';
          document.body.appendChild(notification);
        }

        notification.textContent = message;
        notification.style.backgroundColor = type === 'error' ? '#e74c3c' : '#27ae60';
        notification.style.display = 'block';

        setTimeout(() => {
          notification.style.display = 'none';
        }, 3000);
      }

      function printAlignment() {
        const alignResults = document.querySelector('.alignments-container');
        if (!alignResults) return;

        const printWindow = window.open('', '_blank');
        printWindow.document.write(`
          <html>
            <head>
              <title>Alignement - HGX</title>
              <style>
                body { font-family: 'Courier New', monospace; font-size: 10px; margin: 20px; }
                .alignment-block { margin-bottom: 20px; page-break-inside: avoid; }
                .alignment-title { font-weight: bold; margin-bottom: 10px; }
                .ruler-line { color: #666; }
                pre { white-space: pre-wrap; }
                @media print {
                  body { font-size: 8px; }
                  .alignment-block { page-break-inside: avoid; }
                }
              </style>
            </head>
            <body>
              <h1>Résultats d'alignement - HGX</h1>
              <div>${alignResults.innerHTML}</div>
            </body>
          </html>
        `);
        printWindow.document.close();
        printWindow.print();
      }
    "))
  ),

  tabPanel("Clonages Hors Base",

           # Conteneur options en haut full width
           div(id = "input_container",
               selectInput("carte_xdna", "Choisir une carte .gb :",
                           choices = list.files(xdna_dir, pattern = "\\.gb$", full.names = FALSE)),
               selectInput("seq_files", "Choisir séquence(s) .seq :",
                           choices = list.files(seq_dir, pattern = "\\.seq$", full.names = FALSE),
                           multiple = TRUE),
               actionButton("align_btn", "Lancer l'alignement",
                            style = "background-color: #b22222; color: white; font-weight: bold; border: none; padding: 8px 16px; border-radius: 4px;")
           ),

           # Conteneur résultats full width
           div(id = "results_container",
               verbatimTextOutput("seq_xdna"),
               verbatimTextOutput("seqs_selected"),

               # Boutons d'export
               conditionalPanel(
                 condition = "output.align_results",
                 div(id = "export_buttons",
                     tags$button("📋 Copier", onclick = "copyToClipboard()", class = "export-btn copy-btn"),
                     downloadButton("download_txt", "💾 Télécharger TXT", class = "export-btn")
                 )
               ),

               uiOutput("align_results")
           )
  ),

  tabPanel("TEST",
           h4("Contenu ici...")
  )
))
