# Inclure les scripts
source("global_clonage.R")
source("ui_clonage.R")
source("server_clonage.R")
source("config_clonage.R")

# Lancer l'application
shinyApp(ui = ui_clonage, server = server_clonage)
