library(tidyverse)
library(jsonlite)

params <- fromJSON("params.json")


x_file_metadata <- params$x_file$samples |>
  select(sample_name,group= selcted_group_name )

x_file <- read_tsv(params$x_file$content)|>
  mutate(MS2_name = case_when(is.na(MS2_name)~MS1_name,
                              .default =MS2_name )) |>
  select(MS2_name, x_file_metadata$sample_name) |>
  na.omit() 
x_file <- x_file[!duplicated(x_file$MS2_name),] |>
  column_to_rownames("MS2_name") |>
  t() |> as.data.frame()



y_file_metadata <- params$y_file$samples |>
  select(sample_name,group= selcted_group_name )

y_file <- read_tsv(params$y_file$content)|>
  mutate(MS2_name = case_when(is.na(MS2_name)~MS1_name,
                              .default =MS2_name )) |>
  select(MS2_name, y_file_metadata$sample_name) |>
  na.omit() 
y_file <- y_file[!duplicated(y_file$MS2_name),] |>
  column_to_rownames("MS2_name") |>
  t() |> as.data.frame()

feature_list <- str_split(params$query_feature,",")[[1]]
# feature <- feature_list[1]
feature_list_all <- intersect(colnames(x_file), colnames(y_file))
intersect_feature_list <- intersect(feature_list, feature_list_all)

point_plot <- function(feature){
  # 手动计算 R²
  df <- tibble(
    x = x_file[, feature],
    y = y_file[, feature]
  )
  
  df <- df %>% mutate(
    x = log10(x + 1),  # +1 防止为0
    y = log10(y + 1)
  )
  # cor(df$x, df_log$y, method = "pearson")
  
  fit <- lm(y ~ x, data = df)
  r2_value <- summary(fit)$r.squared
  r2_label <- sprintf("R² = %.3f", r2_value)  # 保留 3 位小数，可改
  
  df %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 2, color = "#2C7BB6", alpha = 0.8) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "#D7191C",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = r2_label,
      hjust = 1.1, vjust = 1.5,
      size = 5,
      fontface = "bold"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      legend.position = "none"
    ) +
    labs(
      x = "X",
      y = "Y",
      title = feature
    )
  ggsave(filename = str_glue("output/{feature}.pdf"))
}


lapply(intersect_feature_list, function(x){
  point_plot(x)
})








