# ==============================================================================
# UI_CLONAGE.R
# ==============================================================================

library(shiny)
library(shinythemes)

ui_clonage <- navbarPage(
  title = div(style = "color: white; font-weight: bold; font-size: 20px;", "HGX"),
  id = "navbar",
  theme = shinytheme("united"),
  windowTitle = "HGX - Analyseur de séquences",

  # ==============================================================================
  # STYLES CSS ET JAVASCRIPT
  # ==============================================================================
  header = tags$head(
    tags$style(HTML("
      /* ===== LAYOUT PRINCIPAL ===== */
      .results-container {
        display: flex;
        gap: 20px;
        align-items: flex-start;
      }

      /* ===== CONTENEUR DE LÉGENDE ===== */
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

      /* ===== CONTENEUR D'ALIGNEMENTS ===== */
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

      /* ===== BLOCS D'ALIGNEMENT ===== */
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

      /* ===== SURLIGNAGE DES MUTATIONS ===== */
      .mutation-highlight {
        background-color: #ffcccc !important;
        border: 1px solid #ff9999 !important;
        border-radius: 2px !important;
      }

      /* ===== INTERFACE DE RECHERCHE ===== */
      .search-container {
        background: #f8f9fa;
        padding: 15px;
        border-radius: 5px;
        border: 1px solid #dee2e6;
        margin-bottom: 10px;
      }

      .search-result {
        background: #e8f5e8;
        padding: 8px;
        border-radius: 4px;
        margin-top: 10px;
        font-family: monospace;
        font-size: 12px;
        border-left: 4px solid #4caf50;
      }

      /* ===== INDICATEUR DE PROGRESSION POUR ALIGNEMENT ===== */
      .alignment-progress {
        background: #e3f2fd;
        border: 2px solid #2196f3;
        border-radius: 8px;
        padding: 20px;
        margin: 20px 0;
        text-align: center;
        animation: pulse-blue 2s infinite;
      }

      @keyframes pulse-blue {
        0% { background-color: #e3f2fd; border-color: #2196f3; }
        50% { background-color: #bbdefb; border-color: #1976d2; }
        100% { background-color: #e3f2fd; border-color: #2196f3; }
      }

      .alignment-spinner {
        display: inline-block;
        width: 20px;
        height: 20px;
        border: 3px solid #f3f3f3;
        border-top: 3px solid #2196f3;
        border-radius: 50%;
        animation: spin 1s linear infinite;
        margin-right: 10px;
      }

      /* ===== RÉDUCTION HAUTEUR ZONE INFORMATIONS ===== */
      .info-section {
        max-height: 250px;
        overflow-y: auto;
        background: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        padding: 10px;
        font-family: 'Courier New', monospace;
        font-size: 10px;
        line-height: 1.2;
        white-space: pre-wrap;
        word-wrap: break-word;
      }

      .info-section::-webkit-scrollbar {
        width: 8px;
      }

      .info-section::-webkit-scrollbar-track {
        background: #f1f1f1;
        border-radius: 4px;
      }

      .info-section::-webkit-scrollbar-thumb {
        background: #c1c1c1;
        border-radius: 4px;
      }

      .info-section::-webkit-scrollbar-thumb:hover {
        background: #a8a8a8;
      }

      /* ===== STYLE POUR BOUTON ALIGNEMENT ACTIF ===== */
      .align-btn-processing {
        background-color: #ff9800 !important;
        animation: pulse-orange 2s infinite;
        pointer-events: none;
      }

      @keyframes pulse-orange {
        0% { background-color: #ff9800; }
        50% { background-color: #f57c00; }
        100% { background-color: #ff9800; }
      }

      /* ===== BOUTONS D'ACTION CÔTE À CÔTE ===== */
      .action-buttons {
        display: flex;
        gap: 10px;
        justify-content: center;
        align-items: center;
        margin: 20px 0;
      }

      /* ===== ANIMATIONS ===== */
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }

      /* ===== TOOLTIPS POUR NUCLÉOTIDES ===== */
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

    # ==============================================================================
    # JAVASCRIPT POUR INTERACTIONS
    # ==============================================================================
    tags$script(HTML("
      // Configuration centralisée
      const buttonConfig = {
        '#search_seq_btn': {
          processing: {
            class: 'btn-searching',
            text: '🔍 Recherche...',
            showStop: true
          },
          reset: {
            class: '',
            text: '🔍 Rechercher',
            showStop: false
          }
        },
        '#align_btn': {
          processing: {
            class: 'align-btn-processing',
            text: '<span class=\"alignment-spinner\"></span>⏳ Alignement en cours...'
          },
          reset: {
            class: '',
            text: '🔬 Lancer l\\'alignement'
          }
        }
      };

      // Fonction générique
      function setButtonState(buttonId, state) {
        const button = $(buttonId);
        button.removeClass().addClass(state.class || '');
        button.html(state.text);
        button.prop('disabled', state.disabled || false);

        if (state.showStop) {
          $('.stop-search-btn').show();
        }
      }

      // Gestionnaires d'événements simplifiés
      $(document).on('click', '#search_seq_btn', function() {
        setButtonState('#search_seq_btn', buttonConfig['#search_seq_btn'].processing);
      });

      $(document).on('click', '#align_btn', function() {
        setButtonState('#align_btn', buttonConfig['#align_btn'].processing);
      });

      // Gestionnaire unifié des réponses Shiny
      $(document).on('shiny:value', function(event) {
        if (event.name === 'search_in_progress' && event.value === false) {
          setButtonState('#search_seq_btn', buttonConfig['#search_seq_btn'].reset);
          $('.stop-search-btn').hide();
        }

        if (event.name === 'alignment_complete' && event.value === true) {
          setTimeout(() => {
            setButtonState('#align_btn', buttonConfig['#align_btn'].reset);
          }, 1000);
        }
      });
    "))
  ),

  # ==============================================================================
  # ONGLET PRINCIPAL - CLONAGES HORS BASE
  # ==============================================================================
  tabPanel("Clonages Hors Base",
           div(style = "width: 100%; margin: 20px auto; padding: 20px;",

               # ==============================================================================
               # SECTION FICHIERS D'ENTRÉE
               # ==============================================================================
               wellPanel(
                 h4("📁 Fichiers d'entrée", style = "color: #b22222; margin-top: 0;"),

                 # Sélection de la carte GenBank de référence
                 fluidRow(
                   column(12,
                          actionButton("refresh_files", "🔄",
                                       style = "float: right; margin-bottom: 5px; background: transparent; border: 1px solid #ccc;",
                                       title = "Rafraîchir la liste des fichiers GenBank"),
                          selectInput("carte_xdna",
                                      "Carte de référence (.gb):",
                                      choices = NULL,
                                      width = "100%")
                   )
                 ),

                 # Interface de recherche de fichiers .seq
                 div(class = "search-container",
                     h5("🔍 Recherche de fichiers .seq", style = "color: #b22222; margin-top: 0; margin-bottom: 15px;"),

                     # Champs de recherche
                     fluidRow(
                       column(6,
                              textInput("plate_keyword",
                                        label = "Nom de plaque (ex: AU83940):",
                                        value = "",
                                        placeholder = "Tapez le nom de la plaque...",
                                        width = "100%")
                       ),
                       column(6,
                              textInput("seq_keyword",
                                        label = "Mot-clé pour fichiers .seq:",
                                        value = "",
                                        placeholder = "Tapez le mot-clé...",
                                        width = "100%")
                       )
                     ),

                     # Boutons d'action
                     fluidRow(
                       column(8,
                              actionButton("search_seq_btn", "🔍 Rechercher",
                                           style = "background-color: #17a2b8; color: white; border: none; padding: 8px 20px; border-radius: 4px; margin-top: 10px;")
                       ),
                       column(4,
                              actionButton("stop_search_btn", "⏹️ Arrêter",
                                           style = "background-color: #dc3545; color: white; border: none; padding: 8px 20px; border-radius: 4px; margin-top: 10px; display: none;",
                                           class = "stop-search-btn")
                       )
                     ),

                     # Affichage des résultats de recherche
                     conditionalPanel(
                       condition = "output.search_results",
                       div(class = "search-result",
                           htmlOutput("search_results")
                       )
                     )
                 ),

                 # Sélection finale des fichiers .seq trouvés
                 conditionalPanel(
                   condition = "output.seq_files_found",
                   # Affichage dynamique des groupes
                   uiOutput("groups_selection_ui"),

                   # Affichage du résumé de sélection seulement
                   div(id = "selection_summary",
                       style = "margin-top: 10px; padding: 8px; background: #e8f4f8; border-radius: 4px; font-family: monospace; font-size: 12px;",
                       textOutput("selection_summary_text")
                   )
                 )
               ),

               # ==============================================================================
               # SECTION SITES DE RESTRICTION
               # ==============================================================================
               wellPanel(
                 h4("🧬 Sites de restriction", style = "color: #b22222; margin-top: 0;"),

                 # Enzyme 1
                 fluidRow(
                   column(6,
                          selectInput("enzyme1",
                                      "Enzyme 1:",
                                      choices = c("Aucune" = "", "Oligo ou autre bout de séquence" = "CUSTOM", names(get_restriction_enzymes())),
                                      selected = "",
                                      width = "100%")
                   ),
                   column(6,
                          # Input conditionnel pour séquence personnalisée
                          conditionalPanel(
                            condition = "input.enzyme1 == 'CUSTOM'",
                            textInput("enzyme1_custom_seq",
                                      "Séquence oligo1:",
                                      value = "",
                                      placeholder = "ex: GAATTC, GCTAGC, ATGCGATCG...",
                                      width = "100%")
                          ),
                          # Nom optionnel pour la séquence personnalisée
                          conditionalPanel(
                            condition = "input.enzyme1 == 'CUSTOM'",
                            textInput("enzyme1_custom_name",
                                      "Nom séquence 1 (optionnel):",
                                      value = "",
                                      placeholder = "ex: Oligo1",
                                      width = "100%")
                          )
                   )
                 ),

                 # Enzyme 2
                 fluidRow(
                   column(6,
                          selectInput("enzyme2",
                                      "Enzyme 2:",
                                      choices = c("Aucune" = "", "Oligo ou autre bout de séquence" = "CUSTOM", names(get_restriction_enzymes())),
                                      selected = "",
                                      width = "100%")
                   ),
                   column(6,
                          # Input conditionnel pour séquence personnalisée
                          conditionalPanel(
                            condition = "input.enzyme2 == 'CUSTOM'",
                            textInput("enzyme2_custom_seq",
                                      "Séquence oligo2:",
                                      value = "",
                                      placeholder = "ex: AAGCTT, CTGCAG, TGCAGTCGA...",
                                      width = "100%")
                          ),
                          # Nom optionnel pour la séquence personnalisée
                          conditionalPanel(
                            condition = "input.enzyme2 == 'CUSTOM'",
                            textInput("enzyme2_custom_name",
                                      "Nom séquence 2 (optionnel):",
                                      value = "",
                                      placeholder = "ex: Oligo2",
                                      width = "100%")
                          )
                   )
                 ),

                 # Affichage des informations sur les sites trouvés
                 div(style = "background: #e8f4f8; padding: 8px; border-radius: 4px; margin-top: 10px; font-family: monospace; font-size: 12px;",
                     textOutput("restriction_info")
                 )

                ),



               # ==============================================================================
               # SECTION FICHIERS AB1
               # ==============================================================================
               conditionalPanel(
                 condition = "output.seq_files_found",
                 wellPanel(
                   h4("📊 Fichiers AB1 correspondants", style = "color: #28a745; margin-top: 0;"),

                   # Zone d'affichage des boutons AB1 individuels (sans bouton de recherche)
                   uiOutput("ab1_buttons_ui")
                 )
               ),

               # ==============================================================================
               # BOUTONS D'ACTION PRINCIPAUX
               # ==============================================================================
               div(class = "action-buttons",
                   actionButton("align_btn", "🔬 Lancer l'alignement",
                                style = "background-color: #b22222; color: white; font-size: 16px; padding: 12px 30px; border: none; border-radius: 5px; font-weight: bold;"),

                   # Bouton d'export HTML conditionnel
                   conditionalPanel(
                     condition = "output.align_results",
                     downloadButton("download_html", "🌐 Télécharger HTML",
                                    style = "background-color: #e74c3c; color: white; border: none; padding: 12px 20px; border-radius: 5px; font-size: 14px;",
                                    title = "Télécharger le rapport complet en HTML")
                   )
               ),

               # ==============================================================================
               # INDICATEUR DE PROGRESSION ALIGNEMENT
               # ==============================================================================
               conditionalPanel(
                 condition = "output.alignment_in_progress",
                 div(class = "alignment-progress",
                     div(class = "alignment-spinner"),
                     tags$strong("Alignement en cours..."),
                     br(),
                     tags$small("Analyse des séquences et génération des résultats colorés")
                 )
               ),

               # ==============================================================================
               # SECTION INFORMATIONS
               # ==============================================================================
               wellPanel(
                 h4("📊 Informations", style = "color: #b22222; margin-top: 0; margin-bottom: 10px;"),

                 # Zone d'affichage compacte avec scroll
                 div(class = "info-section",
                     verbatimTextOutput("seq_xdna_compact"),
                     verbatimTextOutput("seqs_selected_compact")
                 ),

                 # Indicateur de statut
                 div(id = "status-indicator",
                     style = "margin-top: 10px; padding: 5px; text-align: center; font-size: 12px;",
                     uiOutput("processing_status")
                 )
               ),

               # ==============================================================================
               # INSTRUCTIONS UTILISATEUR
               # ==============================================================================
               conditionalPanel(
                 condition = "output.align_results",
                 div(style = "background: #e8f5e8; padding: 10px; border-radius: 4px; margin: 10px 0; border-left: 4px solid #4caf50;",
                     tags$p(style = "margin: 0; font-size: 13px; color: #2e7d32;",
                            "💡 ", tags$strong("Astuce:"), " Survolez les nucléotides dans l'alignement pour voir leur position exacte dans la séquence !")
                 )
               ),

               # ==============================================================================
               # ZONE D'AFFICHAGE DES RÉSULTATS D'ALIGNEMENT
               # ==============================================================================
               uiOutput("align_results")
           )
  ),

  # ==============================================================================
  # ONGLET TEST (POUR DÉVELOPPEMENTS FUTURS)
  # ==============================================================================
  tabPanel("TEST",
           div(style = "padding: 20px;",
               h4("Zone de test", style = "color: #b22222;"),
               p("Cet onglet est réservé pour les développements futurs et les tests.")
           )
  )
)
