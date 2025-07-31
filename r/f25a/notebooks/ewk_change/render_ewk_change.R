library(here)
library(quarto)
library(yaml)

template_type <- "ewk_change"
template_dir <- here::here("notebooks", template_type)
params_path <- file.path(template_dir, paste0(template_type, "_params.yml"))
template_path <- file.path(template_dir, paste0(template_type, "_template.qmd"))

todo <- yaml::read_yaml(params_path)
for (i in seq_along(todo)) {
  params <- todo[[i]]
  quarto::quarto_render(
    input = template_path,
    output_file = paste0(params$fstem, ".pdf"),
    output_format = "pdf",
    execute_params = params,
    quarto_args = c("--output-dir", file.path("_output", template_type))
  )
}
