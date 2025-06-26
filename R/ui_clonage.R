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

      /* Full width results container */
      #results_container {
        width: 100%;
        padding: 15px;
        box-sizing: border-box;
      }

      /* Style bloc alignement */
      #align_results {
        font-family: 'Courier New', monospace !important;
        background: #f8f8f8;
        border: 1px solid #ddd;
        padding: 10px;
        white-space: pre;
        overflow-x: auto;
        max-width: 100%;
        max-height: 400px;
        line-height: 1.2 !important;
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
               actionButton("align_btn", "Lancer l'alignement")
           ),

           # Conteneur résultats full width
           div(id = "results_container",
               verbatimTextOutput("seq_xdna"),
               verbatimTextOutput("seqs_selected"),
               uiOutput("align_results")
           )
  ),

  tabPanel("TEST",
           h4("Contenu ici...")
  )
))
