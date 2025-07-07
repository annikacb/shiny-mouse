FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN R -e '\
install.packages("remotes");\
remotes::install_version("shiny", version = "1.9.1");\
remotes::install_version("shinyWidgets", version = "0.8.6");\
remotes::install_version("DT", version = "0.33");\
remotes::install_version("plotly", version = "4.10.4");\
remotes::install_version("shinydashboard", version = "0.7.2");\
remotes::install_version("ggplot2", version = "3.5.1");\
remotes::install_version("dplyr", version = "1.1.4");\
remotes::install_version("bslib", version = "0.8.0");\
remotes::install_version("shinyjqui", version = "0.4.1");\
remotes::install_github("Schwenk-Lab/ProtPQN");\
remotes::install_version("shinythemes", version = "1.2.0");\
remotes::install_version("tidyr", version = "1.3.1");\
remotes::install_version("heatmaply", version = "1.5.0");\
remotes::install_version("tibble", version = "3.2.1");\
remotes::install_version("stringr", version = "1.5.1");\
remotes::install_version("patchwork", version = "1.2.0");\
'

RUN rm -rf /srv/shiny-server/*
COPY /app/ /srv/shiny-server/

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]