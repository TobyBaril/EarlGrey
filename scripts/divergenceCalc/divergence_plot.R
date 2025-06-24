
library(optparse)

option_list <- list(
  make_option(c("-s", "--species_name"), default=NA, type = "character",
                help="Species name (required)"),
  make_option(c("-g", "--in_gff"), default=NA, type = "character",
                help="GFF with Kimura distances (required)"),
  make_option(c("-o", "--out_directory"), default=NA, type = "character",
                help="Directory to write plots to (required)"),
  make_option(c("-f", "--axis_flip"), action = "store_true",
                default=TRUE, type = "logical", help="Flip x-axis on plots")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check variables are set
if(is.na(opt$species_name)){
  stop("Species name must be supplied")
}
if(is.na(opt$in_gff)){
  stop("Path to input gff must be supplied")
}
if(is.na(opt$out_directory)){
  stop("Path to output directory must be supplied")
}

# Check output folder exists or create
if(!dir.exists(opt$out_directory)){
  message("Creating output directory")
  dir.create(opt$out_directory, recursive = TRUE)
}
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggtext))

# set variables
species_name=opt$species_name
in_gff=opt$in_gff
out_directory=opt$out_directory

# Created plot title
plot_title <- paste0("Repeat landscape of *", gsub("_", " ", species_name), "*")
title_plot <- ggplot() + labs(title = plot_title) + theme(plot.title = element_markdown(hjust = 0.5)) + theme(panel.background = element_blank())

# Read in data
divergence_eg_gff <- read_gff(in_gff)

# Breakdown classification of repeats
divergence_eg_tes_gff <- divergence_eg_gff %>%
  dplyr::mutate(subclass = sub("/.*", "", type),
                superfamily = sub("-.*", "", sub(".*/", "", type)))

# Fix Penelopes
divergence_eg_tes_gff <- divergence_eg_tes_gff %>%
  dplyr::mutate(subclass = ifelse(superfamily == "Penelope", "PLE", subclass)) %>%
  dplyr::mutate(subclass = ifelse(subclass %in% c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Unknown"), subclass, "Other")) %>%
  dplyr::mutate(named_subclass = case_when(subclass == "DNA" ~ "DNA Transposon",
                                           subclass == "LTR" ~ "LTR Retrotransposon",
                                           subclass == "PLE" ~ "Penelope",
                                           subclass == "RC" ~ "Rolling Circle",
                                           .default = subclass))

# Create summary of families
summary_table <- divergence_eg_tes_gff %>%
  as_tibble %>%
  group_by(Name) %>%
  dplyr::select(width, KIMURA80, subclass, superfamily, Name) %>%
  mutate(total_bp = sum(width),
         mean_width = round(mean(width), digits = 2),
         min_width = round(min(width), digits = 2),
         max_width = round(max(width), digits = 2),
         sd_width = round(sd(width), digits = 2),
         div = as.numeric(KIMURA80),
         mean_div = round(mean(div), digits = 2),
         min_div = round(min(div), digits = 2),
         max_div = round(max(div), digits = 2),
         sd_div = round(sd(div), digits = 2)
  ) %>%
  dplyr::select(-width, -div, -KIMURA80) %>%
  base::unique() %>%
  dplyr::arrange(subclass, superfamily, Name) %>%
  dplyr::rename(family = Name)

# Save summary to file
readr::write_tsv(x = summary_table, file = paste0(out_directory, "/", species_name, "_summary_table.tsv"))

# Sum lengths to create data for plots (remove subclasses not in standard set and repeats which Kimura was not calculated for)
divergence_eg_tes_rounded_for_plot  <- divergence_eg_tes_gff %>%
  filter(!is.na(KIMURA80)) %>%
  mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
  dplyr::filter(KIMURA80 <= 0.5) %>%
  as_tibble() %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2)) %>%
  group_by(named_subclass, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup() %>%
  dplyr::select(subclass, named_subclass, KIMURA80, KIMURA_SUM) %>%
  base::unique() %>%
  arrange(named_subclass, KIMURA80)

# Set fill colours
fill_colours <- tibble(subclass = c("DNA", "LINE", "LTR", "PLE", "RC", "SINE", "Other", "Unknown"),
                       named_subclass = c("DNA Transposon", "LINE", "LTR Retrotransposon", "Penelope", "Rolling Circle", "SINE", "Other", "Unknown"),
                       fill_colour = c("#E32017", "#0098D4", "#00782A", "#7156A5", "#EE7C0E", "#9B0056", "#F3A9BB", "#A0A5A9")) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot$subclass) %>%
  arrange(named_subclass) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot$subclass)

# Create and save main plots
kimura_plot <- ggplot(divergence_eg_tes_rounded_for_plot,
                      aes(x = KIMURA80, y = KIMURA_SUM, fill = named_subclass)) +
  geom_col(position = "stack", width = 0.01) +
  theme_bw() +
  labs(title = plot_title) + theme(plot.title = element_markdown(hjust = 0.5)) +
  scale_fill_manual(values = fill_colours$fill_colour, name = "TE Subclass")
subclass_kimura_plot <- kimura_plot + scale_y_continuous(expand = c(0.01,0), name = "Base pairs")
# Flip axis if desired
if(opt$axis_flip == TRUE){
  subclass_kimura_plot <- subclass_kimura_plot +
    scale_x_continuous(limits = c(0.51, -0.01),
                       expand = c(0,0), name = "Kimura 2-Parameter Distance",
                       trans = "reverse")
} else {
  subclass_kimura_plot <- subclass_kimura_plot +
    scale_x_continuous(limits = c(-0.01, 0.51),
                       expand = c(0,0), name = "Kimura 2-Parameter Distance")
}
ggsave(plot = subclass_kimura_plot, filename = paste0(out_directory, "/", species_name, "_classification_landscape.pdf"), device = "pdf", width = 12.85, height = 8.5)

split_subclass_kimura_plot <- kimura_plot + scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) + facet_grid(subclass~., scales = "free")
# Flip axis if desired
if(opt$axis_flip == TRUE){
  split_subclass_kimura_plot <- split_subclass_kimura_plot +
    scale_x_continuous(limits = c(0.51, -0.01),
                       expand = c(0,0), name = "Kimura 2-Parameter Distance",
                       trans = "reverse")
} else {
  split_subclass_kimura_plot <- split_subclass_kimura_plot +
    scale_x_continuous(limits = c(-0.01, 0.51),
                       expand = c(0,0), name = "Kimura 2-Parameter Distance")
}
ggsave(plot = split_subclass_kimura_plot, filename = paste0(out_directory, "/", species_name, "_split_class_landscape.pdf"), device = "pdf", width = 12.85, height = 8.5)

# Perform maths for more divided plot
divergence_eg_tes_rounded_for_superfamily_plot  <- divergence_eg_tes_gff %>%
  as_tibble() %>%
  filter(!is.na(KIMURA80)) %>%
  mutate(KIMURA80 = as.numeric(KIMURA80)) %>%
  dplyr::filter(KIMURA80 <= 0.5) %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2),
                type = sub("-.*", "", type)) %>%
  group_by(superfamily, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup() %>%
  dplyr::select(type, subclass, superfamily, KIMURA80, KIMURA_SUM) %>%
  base::unique() %>%
  arrange(subclass, superfamily, KIMURA80)

# Split data as necessary
divergence_eg_tes_rounded_for_superfamily_plot <- split(divergence_eg_tes_rounded_for_superfamily_plot,
                                                        f = divergence_eg_tes_rounded_for_superfamily_plot$subclass)

# Create plots of superfamilies of DNA transposons, LINEs, LTR retrotransposons and SINEs
kimura_superfamily_plot_1 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$DNA,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3))
if (inherits(try(ggplot_build(kimura_superfamily_plot_1)), "try-error")) 
  kimura_superfamily_plot_1 <- NULL

kimura_superfamily_plot_2 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$LINE,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_brewer(palette = "Blues", direction = -1)
if (inherits(try(ggplot_build(kimura_superfamily_plot_2)), "try-error")) 
  kimura_superfamily_plot_2 <- NULL
                     
kimura_superfamily_plot_3 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$LTR,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_brewer(palette = "Greens", direction = -1)
if (inherits(try(ggplot_build(kimura_superfamily_plot_3)), "try-error")) 
  kimura_superfamily_plot_3 <- NULL

kimura_superfamily_plot_4 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$SINE,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1)
if (inherits(try(ggplot_build(kimura_superfamily_plot_4)), "try-error")) 
  kimura_superfamily_plot_4 <- NULL

# flip axis if desired
# Flip axis if desired
if(opt$axis_flip == TRUE){
  kimura_superfamily_plot_1 <- kimura_superfamily_plot_1 +
    scale_x_continuous(limits = c(0.51, -0.01),
                       expand = c(0,0), name = "", trans = "reverse")
  kimura_superfamily_plot_2 <- kimura_superfamily_plot_2 +
    scale_x_continuous(limits = c(0.51, -0.01),
                       expand = c(0,0), name = "", trans = "reverse")
  kimura_superfamily_plot_3 <- kimura_superfamily_plot_3 +
    scale_x_continuous(limits = c(0.51, -0.01),
                       expand = c(0,0), name = "", trans = "reverse")
  kimura_superfamily_plot_4 <- kimura_superfamily_plot_4 +
    scale_x_continuous(limits = c(0.51, -0.01), trans = "reverse",
                       expand = c(0,0), name = "Kimura 2-Parameter Distance")
} else {
  kimura_superfamily_plot_1 <- kimura_superfamily_plot_1 +
    scale_x_continuous(limits = c(-0.01, 0.51),
                       expand = c(0,0), name = "")
  kimura_superfamily_plot_2 <- kimura_superfamily_plot_2 +
    scale_x_continuous(limits = c(-0.01, 0.51),
                       expand = c(0,0), name = "")
  kimura_superfamily_plot_3 <- kimura_superfamily_plot_3 +
    scale_x_continuous(limits = c(-0.01, 0.51),
                       expand = c(0,0), name = "")
  kimura_superfamily_plot_4 <- kimura_superfamily_plot_4 +
    scale_x_continuous(limits = c(-0.01, 0.51),
                       expand = c(0,0), name = "Kimura 2-Parameter Distance")
}

# Combine plots and title
superfamily_kimura_plot <- plot_grid(kimura_superfamily_plot_1, kimura_superfamily_plot_2, kimura_superfamily_plot_3, kimura_superfamily_plot_4, 
                                     ncol = 1, align = "v")
superfamily_kimura_plot_titled <- plot_grid(title_plot, superfamily_kimura_plot, ncol = 1, rel_heights = c(1, 30))

# Save divided plot
ggsave(plot = superfamily_kimura_plot_titled, filename = paste0(out_directory, "/", species_name, "_superfamily_div_plot.pdf"), device = "pdf", width = 12.85, height = 8.5)
