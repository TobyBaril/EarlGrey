library(optparse)

option_list <- list(
  make_option(c("-s", "--species_name"), default=NA, type = "character", help="Species name (required)"),
  make_option(c("-g", "--in_gff"), default=NA, type = "character", help="GFF with Kimura distances (required)"),
  make_option(c("-o", "--out_directory"), default=NA, type = "character", help="Directory to write plots to (required)")
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

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggtext))

# Created plot title
plot_title <- paste0("Repeat landscape of *", gsub("_", " ", opt$species_name), "*")
title_plot <- ggplot() + labs(title = plot_title) + theme(plot.title = element_markdown(hjust = 0.5)) + theme(panel.background = element_blank())

# Read in data, remove repeats which Kimura was not calculated for
divergence_eg_gff <- read_gff(opt$in_gff) %>%
  mutate(KIMURA80 = as.double(KIMURA80)) %>%
  filter(!is.na(KIMURA80), KIMURA80 <= 0.5)

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
                                           subclass == "Unknown" ~ "Unclassified",
                                           .default = subclass))

# Sum lengths to create data for plots (remove subclasses not in standard set)
divergence_eg_tes_rounded_for_plot  <- divergence_eg_tes_gff %>%
  as_tibble() %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2)) %>%
  group_by(named_subclass, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup() %>%
  dplyr::select(subclass, named_subclass, KIMURA80, KIMURA_SUM) %>%
  base::unique() %>%
  arrange(named_subclass, KIMURA80)

divergence_eg_tes_rounded_for_plot$named_subclass %<>% 
  as.factor() %>%
  ordered(levels = c("DNA Transposon", "Rolling Circle", "Penelope", "LINE", "SINE", "LTR Retrotransposon", "Other (Simple Repeat, Microsatellite, RNA)", "Unclassified"))

# Set fill colours
fill_colours <- data.frame(subclass = c("DNA", "RC", "PLE", "LINE", "SINE", "LTR", "Other", "Unknown"),
                       named_subclass = c("DNA Transposon", "Rolling Circle", "Penelope", "LINE", "SINE", "LTR Retrotransposon", "Other (Simple Repeat, Microsatellite, RNA)", "Unclassified"),
                       fill_colour = c("#E32017", "#EE7C0E", "#7156A5", "#0098D4", "#9B0056", "#00782A", "#F3A9BB", "#A0A5A9")) %>%
  filter(subclass %in% divergence_eg_tes_rounded_for_plot$subclass) 

col <- fill_colours$fill_colour
names(col) <- fill_colours$named_subclass

# Create and save main plots
kimura_plot <- ggplot(divergence_eg_tes_rounded_for_plot,
                      aes(x = KIMURA80, y = KIMURA_SUM, fill = named_subclass)) +
  geom_col(position = "stack", width = 0.01) +
  scale_x_reverse(limits = c(1,0), expand = c(0,0), name = "Kimura 2-Parameter Distance") +
  theme_classic() +
  labs(title = plot_title) + 
  theme(plot.title = element_markdown(hjust = 0.5)) +
  scale_fill_manual(values = col, name = "TE Classification")

subclass_kimura_plot <- kimura_plot + scale_y_continuous(expand = c(0.01,0), name = "Base pairs")

ggsave(plot = subclass_kimura_plot, 
       filename = paste0(opt$out_directory, "/", opt$species_name, "_classification_landscape.pdf"),
       device = "pdf", 
       scale = 1,
       width = 297, 
       height = 210,
       units = "mm",
       dpi = 300,
       limitsize = FALSE)

split_subclass_kimura_plot <- kimura_plot + scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) + facet_grid(subclass~., scales = "free")

ggsave(plot = split_subclass_kimura_plot, 
       filename = paste0(opt$out_directory, "/", opt$species_name, "_split_class_landscape.pdf"), 
       device = "pdf", 
       scale = 1,
       width = 297, 
       height = 210,
       units = "mm",
       dpi = 300,
       limitsize = FALSE)

# Perform maths for more divided plot
divergence_eg_tes_rounded_for_superfamily_plot  <- divergence_eg_tes_gff %>%
  as_tibble() %>%
  dplyr::mutate(KIMURA80 = round(x = KIMURA80, digits = 2),
                type = sub("-.*", "", type)) %>%
  group_by(superfamily, KIMURA80) %>%
  mutate(KIMURA_SUM = sum(width)) %>%
  ungroup()%>%
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
  scale_x_reverse(limits = c(1,0), expand = c(0,0), name = "") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3))
if (inherits(try(ggplot_build(kimura_superfamily_plot_1)), "try-error")) 
  kimura_superfamily_plot_1 <- ggplot()

kimura_superfamily_plot_2 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$LINE,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  scale_x_reverse(limits = c(1,0), expand = c(0,0), name = "") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_brewer(palette = "Blues", direction = -1)
if (inherits(try(ggplot_build(kimura_superfamily_plot_2)), "try-error")) 
  kimura_superfamily_plot_2 <- ggplot()

kimura_superfamily_plot_3 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$LTR,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  scale_x_reverse(limits = c(1,0), expand = c(0,0), name = "") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_brewer(palette = "Greens", direction = -1)
if (inherits(try(ggplot_build(kimura_superfamily_plot_3)), "try-error")) 
  kimura_superfamily_plot_3 <- ggplot()

kimura_superfamily_plot_4 <- ggplot(divergence_eg_tes_rounded_for_superfamily_plot$SINE,
                                    aes(x = KIMURA80, y = KIMURA_SUM, fill = superfamily)) +
  geom_col(position = "stack", width = 0.01, colour = "black", linewidth = 0.2) +
  scale_x_reverse(limits = c(1,0), expand = c(0,0), name = "Kimura 2-Parameter Distance") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  scale_y_continuous(name = "Base pairs", labels = function(x) format(x, scientific = TRUE)) +
  facet_grid(subclass~., scales = "free") +
  guides(fill=guide_legend(ncol=3)) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1)
if (inherits(try(ggplot_build(kimura_superfamily_plot_4)), "try-error")) 
  kimura_superfamily_plot_4 <- ggplot()

# Combine plots and title
superfamily_kimura_plot <- plot_grid(kimura_superfamily_plot_1, kimura_superfamily_plot_2, kimura_superfamily_plot_3, kimura_superfamily_plot_4, 
          ncol = 1, align = "v")
superfamily_kimura_plot_titled <- plot_grid(title_plot, superfamily_kimura_plot, ncol = 1, rel_heights = c(1, 30))

# Save divided plot
ggsave(plot = superfamily_kimura_plot_titled, 
       filename = paste0(opt$out_directory, "/", opt$species_name, "_superfamily_div_plot.pdf"), 
       device = "pdf", 
       scale = 1,
       width = 297, 
       height = 210,
       units = "mm",
       dpi = 300,
       limitsize = FALSE)
