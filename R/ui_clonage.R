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

      /* Style bloc alignement avec scroll horizontal */
      #align_results {
        font-family: 'Courier New', monospace !important;
        background: #f8f8f8;
        border: 1px solid #ddd;
        padding: 15px;
        white-space: pre;
        overflow-x: auto;
        overflow-y: auto;
        max-width: 100%;
        max-height: 600px;
        line-height: 1.4 !important;
        border-radius: 5px;
        position: relative;
      }

      /* Style pour la légende des couleurs - FIXE en haut */
      #legend-fixed {
        position: sticky !important;
        top: 0 !important;
        z-index: 1000 !important;
        margin-bottom: 15px !important;
        padding: 10px !important;
        background: #f0f0f0 !important;
        border: 1px solid #ddd !important;
        border-radius: 5px !important;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1) !important;
      }

      .color-legend {
        margin-bottom: 20px;
        padding: 15px;
        background: #e8f4fd;
        border: 1px solid #bee5eb;
        border-radius: 8px;
        font-family: Arial, sans-serif;
      }

      .color-legend b {
        color: #0c5460;
        font-size: 14px;
      }

      /* Style pour les alignements individuels */
      .alignment-block {
        margin-bottom: 25px;
        padding: 15px;
        background: #ffffff;
        border: 1px solid #e0e0e0;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }

      .alignment-title {
        color: #2c3e50;
        font-weight: bold;
        margin-bottom: 10px;
        font-family: Arial, sans-serif;
        border-bottom: 2px solid #3498db;
        padding-bottom: 5px;
      }

      .alignment-score {
        color: #7f8c8d;
        font-style: italic;
        margin-top: 8px;
        font-family: Arial, sans-serif;
      }

      /* Style pour la règle de numérotation */
      .ruler-line {
        color: #666 !important;
        font-weight: normal !important;
        border-bottom: 1px solid #ddd;
        margin-bottom: 2px;
        padding-bottom: 2px;
        background: #f9f9f9;
      }

      /* Amélioration des tooltips */
      span[title]:hover {
        background-color: #ffffcc;
        border-radius: 2px;
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

      /* Assurer que la police monospace est bien appliquée partout */
      #align_results * {
        font-family: 'Courier New', monospace !important;
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
    ")),

    # JavaScript pour la fonctionnalité de copie
    tags$script(HTML("
      function copyToClipboard() {
        const alignResults = document.getElementById('align_results');
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
        const alignResults = document.getElementById('align_results');
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
