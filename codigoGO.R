library(tidyverse)
library(stringr)

home_dir <- Sys.getenv("HOME")

path_go <- file.path(
  home_dir,
  "proyecto",
  "trabajoR",
  "TerminosGO"
)




leer_GO <- function(archivo){
  
  read.delim(
    file.path(path_go, archivo),
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    
    dplyr::select(
      GO = Name,
      Genes = `Result count`,
      Benjamini
    ) %>%
    
    mutate(
      Genes = as.numeric(Genes),
      Benjamini = as.numeric(Benjamini)
    )
  
}


bp_up   <- leer_GO("SobreExpresados Biologico.tsv")
bp_down <- leer_GO("InfraExpresados Biologico.tsv")

mf_up   <- leer_GO("SobreExpresados Molecular.tsv")
mf_down <- leer_GO("InfraExpresados Molecular.tsv")

cc_up   <- leer_GO("SobreExpresados Celular.tsv")
cc_down <- leer_GO("InfraExpresados Celular.tsv")




combinar_GO <- function(up, down, categoria){
  
  full_join(
    up,
    down,
    by = "GO",
    suffix = c("_Up", "_Down")
  ) %>%
    
    replace_na(
      list(
        Genes_Up = 0,
        Genes_Down = 0,
        Benjamini_Up = 1,
        Benjamini_Down = 1
      )
    ) %>%
    
    mutate(
      
      Categoria = categoria,
      
      Total = Genes_Up + Genes_Down
      
    )
  
}



bp <- combinar_GO(
  bp_up,
  bp_down,
  "Biological Process"
)

mf <- combinar_GO(
  mf_up,
  mf_down,
  "Molecular Function"
)

cc <- combinar_GO(
  cc_up,
  cc_down,
  "Cellular Component"
)

datos <- bind_rows(bp, mf, cc)





datos <- datos %>%
  
  group_by(Categoria) %>%
  
  slice_max(
    Total,
    n = 12,
    with_ties = FALSE
  ) %>%
  
  ungroup()





datos_plot <- datos %>%
  
  pivot_longer(
    
    cols = c(
      Genes_Up,
      Genes_Down
    ),
    
    names_to = "Regulacion",
    
    values_to = "Genes"
    
  ) %>%
  
  mutate(
    
    Regulacion = recode(
      Regulacion,
      Genes_Up = "Up",
      Genes_Down = "Down"
    ),
    
    Genes = ifelse(
      Regulacion == "Down",
      -Genes,
      Genes
    ),
    
    GO = str_wrap(GO, width = 40)
    
  )







datos_plot <- datos_plot %>%
  
  group_by(Categoria, GO) %>%
  
  mutate(Total = max(abs(Genes))) %>%
  
  ungroup()







datos_plot <- datos %>%
  pivot_longer(
    cols = c(Genes_Up, Genes_Down),
    names_to = "Regulacion",
    values_to = "Genes"
  ) %>%
  mutate(
    Regulacion = recode(Regulacion,
                        Genes_Up = "Up",
                        Genes_Down = "Down"),
    Genes = ifelse(Regulacion == "Down", -Genes, Genes),
    GO = stringr::str_wrap(GO, 35)
  ) %>%
  group_by(Categoria) %>%
  mutate(
    GO = forcats::fct_reorder(GO, Total)
  ) %>%
  ungroup()

p <- ggplot(datos_plot,
            aes(x = Genes,
                y = GO,
                fill = Regulacion)) +

  geom_col(width = 0.75) +

  geom_vline(xintercept = 0,
             linewidth = 0.4,
             colour = "black") +

  facet_grid(Categoria ~ .,
             scales = "free_y",
             space = "free_y") +

  scale_fill_manual(
    values = c(
      Up = "#E64B35",
      Down = "#4DBBD5"
    )
  ) +
  
  scale_x_continuous(
    breaks = seq(-200, 350, by = 50),
    labels = abs,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    x = "Number of genes",
    y = NULL
  ) +

  theme_classic(base_size = 14) +

  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(face = "bold"),
    panel.spacing = unit(0.8, "cm"),
    plot.margin = margin(10, 10, 10, 20)
  )








ggsave(
  filename = file.path(path_go, "GO_barplot_publication.png"),
  plot = p,
  width = 16,
  height = 12,
  dpi = 600
)

ggsave(
  filename = file.path(path_go, "GO_barplot_publication.pdf"),
  plot = p,
  width = 16,
  height = 12
)

print(p)