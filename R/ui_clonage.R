library(shiny)
library(shinythemes)

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
      @media (max-width: 768px) {
        .container-fluid { padding-left: 10px; padding-right: 10px; }
        .sidebarPanel, .mainPanel { width: 100% !important; float: none !important; }
      }
      /* Pour réduire espacement dans align_results */
      #align_results {
        font-family: 'Courier New', monospace;
        line-height: 1.1 !important;
        margin: 0 !important;
        padding: 0 !important;
        white-space: pre-wrap;
      }
    "))
  ),

  tabPanel("Clonages Hors Base",
           sidebarLayout(
             sidebarPanel(
               selectInput("carte_xdna", "Choisir une carte .gb :",
                           choices = list.files(xdna_dir, pattern = "\\.gb$", full.names = FALSE)),
               selectInput("seq_files", "Choisir séquence(s) .seq :",
                           choices = list.files(seq_dir, pattern = "\\.seq$", full.names = FALSE),
                           multiple = TRUE),
               actionButton("align_btn", "Lancer l'alignement")
             ),

             mainPanel(
               verbatimTextOutput("seq_xdna"),
               verbatimTextOutput("seqs_selected"),
               uiOutput("align_results")
             )
           )
  ),

  tabPanel("TEST",
           h4("Contenu ici...")
  )
))
