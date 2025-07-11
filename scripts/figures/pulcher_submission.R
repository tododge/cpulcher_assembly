##################
#### Figure 2 ####
##################

library(tidyverse)
library(pafr)
library(circlize)
library(cowplot)


### Load & Process Alignment
ali <- read_paf("~/pulcher/alignment/cryege_v2024.1_rc.fa_2_crypul_v2024.1.fa.paf")
prim_alignment <- filter(filter_secondary_alignments(ali), alen > 1e6)

rename_chr <- function(x, prefix) gsub("chr", paste0(prefix, "-chr"), x)

t_bed <- prim_alignment %>%
  select(tname, tstart, tend) %>%
  mutate(tname = rename_chr(tname, "crypul"))

q_bed <- prim_alignment %>%
  select(qname, qstart, qend, strand) %>%
  mutate(qname = rename_chr(qname, "cryege"))

chrominfo <- bind_rows(
  prim_alignment %>%
    select(name = tname, len = tlen) %>%
    distinct() %>%
    mutate(name = rename_chr(name, "crypul")),
  prim_alignment %>%
    select(name = qname, len = qlen) %>%
    distinct() %>%
    mutate(name = rename_chr(name, "cryege")) %>%
    arrange(desc(name))
)
colnames(chrominfo) <- c("V1", "V2")

### Load Gaps and Telomeres

load_bed_with_vals <- function(path, prefix, value) {
  df <- read.table(path, sep = "\t", header = FALSE)
  df$vals <- value
  df <- df[grep("chr", df$V1), ]
  df$V1 <- rename_chr(df$V1, prefix)
  df
}

breaks <- bind_rows(
  load_bed_with_vals("~/pulcher/alignment/crypul_v2024.1.fa_scaffold_breaks.txt", "crypul", 1),
  load_bed_with_vals("~/pulcher/alignment/cryege_v2024.1_rc.fa_scaffold_breaks.txt", "cryege", 1)
)

telo <- bind_rows(
  load_bed_with_vals("~/pulcher/alignment/crypul_v2024.1.fa_telo.txt", "crypul", 0.5),
  load_bed_with_vals("~/pulcher/alignment/cryege_v2024.1_rc.fa_telo.txt", "cryege", 0.5)
) %>%
  rename(chr = V1, start = V2, end = V3, value1 = vals) %>%
  select(chr, start, end, value1)

genomeinfo <- read.table("~/pulcher/alignment/crypto_assemblies_genomeinfo.bed", sep = "\t", header = TRUE)
chromosome.index <- c(paste0("cryege-chr", 1:15), rev(paste0("crypul-chr", 1:15)))
genomeinfo[, 1] = factor(genomeinfo[ ,1], levels = chromosome.index)


### Circos Plot

pdf("~/pulcher/figs/circos_plot.pdf", width = 4.5, height = 4.5)
circos.clear()
rcols <- scales::alpha(ifelse(q_bed$strand == "+", "#e3be6b", "#3cb4a4"), alpha=0.8)
circos.par("start.degree" = 87.5, gap.after = c(rep(3, 14), 5, rep(3, 14), 5))

circos.genomicInitialize(genomeinfo, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
              gsub("crypul-chr", "", gsub("cryege-chr", "", CELL_META$sector.index)), cex = 0.6, niceFacing = TRUE)
}, track.height = mm_h(1), cell.padding = c(0, 0, 0, 0), bg.border = NA)

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
}, bg.col="grey90", bg.border=F, track.height = mm_h(2), cell.padding = c(0, 0, 0, 0))

circos.genomicTrack(telo,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 0.6, ...)
                    }, track.index = 2,
                    ylim = c(0,1)
)


circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top", direction = "outside")
})

circos.genomicTrack(breaks,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 1, ybottom = 0)
                    }, track.index = 2,
                    ylim = c(0,1)
)

highlight.chromosome(paste0("crypul-chr", (1:15)), col = "#A69D97", track.index = 1)
highlight.chromosome(paste0("cryege-chr", (1:15)), col = "#69b9cd", track.index = 1)

circos.genomicLink(q_bed, t_bed, col=rcols, border = NA)
dev.off()
circos.clear()


### Dotplot Functions
plot_paf_dotplot <- function(file, seqX, Xcoords, seqY, Ycoords, align_name, offsetX = 0, xbreaks = 10e6, ybreaks = 10e6) {
  mum <- read_paf(file) %>%
    transmute(
      start1 = as.numeric(tstart),
      end1 = as.numeric(tend),
      name1 = tname,
      name2 = qname,
      start2 = as.numeric(ifelse(strand == "+", qstart, qend)),
      end2 = as.numeric(ifelse(strand == "+", qend, qstart))
    )
  
  ggplot(filter(mum, name1 == seqX & name2 == seqY),
         aes(x = start1 + offsetX, y = start2, xend = end1 + offsetX, yend = end2)) +
    geom_segment(lineend = "butt", lwd = 1) +
    theme_bw(base_size = 8) +
    theme(panel.grid = element_blank()) +
    annotation_custom(grid::textGrob(label = align_name, x = 0.1, y = 0.9,
                                     hjust = 0, gp = grid::gpar(fontsize = 8))) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, max(Xcoords), by = xbreaks), labels = ~ . / 1e6) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, max(Ycoords), by = ybreaks), labels = ~ . / 1e6) +
    coord_fixed(xlim = Xcoords, ylim = Ycoords)
}

chrYalign <- plot_paf_dotplot("~/pulcher/alignment/crypul_v2024.1_chrY.fa_2_cryege_v2024.1_ctg_008_Y_32451619_45003156.fa.paf",
                              "ctg_008_Y_32451619_45003156", 32451619 + c(0, 8171628),
                              "chrY", c(0, 6659363),
                              "Y chromosome (partial)", offsetX = 32451619, xbreaks = 1e6, ybreaks = 1e6) +
  xlab(expression(paste(italic("C. egeriae"), " chrY position (Mb)"))) +
  ylab(expression(paste(italic("C. pulcher"), " chrY position (Mb)"))) +
  annotate("rect", xmin = 35e6, xmax = 40623247, ymin = 0, ymax = 0.18e6, fill = "#69b9cd", color = "#092d44", size = 0.25) +
  annotate("rect", xmin = 32451619, xmax = 32451619 + 0.18e6, ymin = 0, ymax = 6659363, fill = "#A69D97", color = "#272e34", size = 0.25)

chrXalign <- plot_paf_dotplot("~/pulcher/alignment/crypul_v2024.1_chr7.fa_2_cryege_v2024.1_chr7.fa.paf",
                              "chr7", c(0, 71111990),
                              "chr7", c(0, 77142432),
                              "X chromosome") +
  xlab(expression(paste(italic("C. egeriae"), " chr7 position (Mb)"))) +
  ylab(expression(paste(italic("C. pulcher"), " chr7 position (Mb)"))) +
  annotate("rect", xmin = 35e6, xmax = 45e6, ymin = 0, ymax = 2e6, fill = "#69b9cd", color = "#092d44", size = 0.25) +
  annotate("rect", xmin = 0, xmax = 2e6, ymin = 37e6, ymax = 48e6, fill = "#A69D97", color = "#272e34", size = 0.25)


### Combine Plots

blank <- ggplot() + theme_void()
sexchrom <- plot_grid(chrXalign, chrYalign, nrow = 2, labels = c("(b)", "(c)"), align = "v", rel_heights = c(1,1), label_size = 10)
circos_abc <- plot_grid(blank, sexchrom, ncol = 2, labels = c("(a)", ""), rel_widths = c(2.2,1), label_size = 10)

ggsave("~/pulcher/figs/figure_circos.pdf", circos_abc, height = 115, width = 169, units = "mm")


##################
#### Figure 3 ####
##################

library(tidyverse)
library(cowplot)

# Helper function: preprocess heterozygosity data
process_het <- function(het_df, chr_len_df) {
  # Join chromosome lengths
  het_df <- left_join(het_df, chr_len_df, by = c("chromosome" = "V1"))
  colnames(het_df)[ncol(het_df)] <- "chr_length"
  
  # Extract numeric chromosome name
  het_df$name <- as.numeric(gsub("\\D", "", het_df$chromosome))
  
  # Compute chromosome offsets and cumulative positions
  chr_offsets <- het_df %>%
    group_by(name) %>%
    summarise(chr_length = max(chr_length), .groups = "drop") %>%
    arrange(-chr_length) %>%
    mutate(tot = cumsum(chr_length) - chr_length)
  
  het_updt <- left_join(chr_offsets, het_df, by = c("name", "chr_length")) %>%
    arrange(name, window_pos_1) %>%
    mutate(
      window_pos_1_cum = window_pos_1 + tot,
      window_pos_2_cum = window_pos_2 + tot,
      name = fct_reorder(as.factor(name), chr_length)
    )
  
  # Compute center positions for axis labels
  axisdf <- het_updt %>%
    group_by(name) %>%
    summarize(
      center = (as.numeric(max(window_pos_2_cum)) + as.numeric(min(window_pos_2_cum))) / 2,
      .groups = "drop"
    ) %>%
    arrange(-center)
  
  list(het_updt = het_updt, axisdf = axisdf)
}

plot_het_rect <- function(het_updt, axisdf, label, fill_colors, species_label) {
  ggplot() +
    geom_rect(data = het_updt, aes(
      xmin = window_pos_1_cum,
      xmax = window_pos_2_cum,
      ymin = 0,
      ymax = avg_pi * 1000,
      fill = name
    )) +
    scale_fill_manual(values = rep(fill_colors, 100)) +
    scale_x_continuous(label = axisdf$name, breaks = axisdf$center) +
    scale_y_continuous(limits = c(0, 20), breaks = seq(0, 15, 5)) +
    theme_classic(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    xlab(species_label) +
    ylab("het. per kb") +
    annotate(
      "text", x = 1e6, y = 18,
      label = label,
      size = 3, hjust = 0, vjust = 0
    )
}

plot_het_hist <- function(het_updt, fill_colors) {
  ggplot(data = het_updt, aes(x = avg_pi * 1000)) +
    geom_histogram(color = fill_colors[1], fill = fill_colors[2], bins = 64, size = 0.25) +
    scale_y_continuous(limits = c(0, 200), breaks = seq(0, 150, 50)) +
    scale_x_continuous(limits = c(-1, 15)) +
    theme_classic(base_size = 8) +
    theme(panel.grid = element_blank()) +
    xlab("het. per kb") +
    ylab("no. bins") +
    annotate(
      "text", x = 15, y = 170,
      label = paste("median=", signif(median(het_updt$avg_pi) * 1000, 3)),
      size = 2.5, hjust = 1
    )
}


### --- Data 1: C. pulcher HiFi ---
chr_len <- read.table("~/pulcher/het/crypul_v2024.1_rmY.fa.fai")[, c(1, 2)]
het_hifi <- read.table("~/pulcher/het/crypul_v2024.1_rmY.fa.hifireads.1000000_pi.txt", header = TRUE)

d1 <- process_het(het_hifi, chr_len)
fig_hifi <- plot_grid(
  plot_het_rect(d1$het_updt, d1$axisdf, "C. pulcher individual #1 heterozygosity (HiFi reads)",
                c("#272e34", "#A69D97"), expression(paste(italic("C. pulcher"), " chromosome"))),
  plot_het_hist(d1$het_updt, c("#272e34", "#A69D97")),
  ncol = 2, labels = c("(a)", "(b)"), rel_widths = c(3, 1), label_size = 10
)

### --- Data 2: C. pulcher Omni-C ---
het_omni <- read.table("~/pulcher/het/crypul_v2024.1_rmY.fa.omniCreads.1000000_pi.txt", header = TRUE)

d2 <- process_het(het_omni, chr_len)
fig_omni <- plot_grid(
  plot_het_rect(d2$het_updt, d2$axisdf, "C. pulcher individual #2 heterozygosity (Omni-C reads)",
                c("#272e34", "#A69D97"), expression(paste(italic("C. pulcher"), " chromosome"))),
  plot_het_hist(d2$het_updt, c("#272e34", "#A69D97")),
  ncol = 2, labels = c("(c)", "(d)"), rel_widths = c(3, 1), label_size = 10
)

### --- Data 3: C. egeriae HiFi ---
chr_len <- read.table("~/pulcher/het/cryege_v2024.1.fa.fai")[, c(1, 2)]
het_hifi <- read.table("~/pulcher/het/cryege_v2024.1_rmY.fa.hifireads.1000000_pi.txt", header = TRUE)

d3 <- process_het(het_hifi, chr_len)
fig_egeriae <- plot_grid(
  plot_het_rect(d3$het_updt, d3$axisdf, "C. egeriae individual #1 heterozygosity (HiFi reads)",
                c("#092d44", "#69b9cd"), expression(paste(italic("C. egeriae"), " chromosome"))),
  plot_het_hist(d3$het_updt, c("#092d44", "#69b9cd")),
  ncol = 2, labels = c("(e)", "(f)"), rel_widths = c(3, 1), label_size = 10
)


####PSMC####

PSMC_all <- read.csv("~/pulcher/psmc/pulcher_egeriae_comb3.txt", header=FALSE, sep = "\t")

PSMC_all <- PSMC_all %>%
  mutate(
    V8 = as.factor(gsub("results3/bootstrap/", "", V8)),
    V6 = trimws(V6),
    V7 = trimws(V7)
  )

# Remove first 6 rows per group
PSMC_filtered <- PSMC_all %>%
  group_by(V6, V7, V8) %>%
  arrange(V1, .by_group = TRUE) %>%
  slice(-1:-6) %>%
  ungroup()

# Define plotting style
style_map <- tribble(
  ~V6,                    ~V7,                ~color,      ~alpha, ~linetype,
  "C. egeriae bootstrap", NA,                 "#69b9cd",    0.5,    "solid",
  "C. egeriae",           NA,                 "#092d44",    1.0,    "solid",
  "C. pulcher bootstrap", "hifi autosomes",   "#A69D97",    0.5,    "solid",
  "C. pulcher",           "hifi autosomes",   "#272e34",    1.0,    "solid",
  "C. pulcher bootstrap", "omniC autosomes",  "#A69D97",    0.5,    "dotted",
  "C. pulcher",           "omniC autosomes",  "#272e34",    1.0,    "dotted"
)

# Start plot
PSMC_plot <- ggplot()

# Add a geom_step layer for each group in style_map
for (i in 1:nrow(style_map)) {
  subset_data <- PSMC_filtered %>%
    filter(
      V6 == style_map$V6[i],
      if (!is.na(style_map$V7[i])) V7 == style_map$V7[i] else TRUE
    )
  
  PSMC_plot <- PSMC_plot +
    geom_step(
      data = subset_data,
      aes(x = V1, y = V2 * 10, group = V8),
      color = style_map$color[i],
      alpha = style_map$alpha[i],
      linetype = style_map$linetype[i]
    )
}

# Finalize styling
PSMC_plot <- PSMC_plot +
  scale_x_log10(
    breaks = c(1e4, 1e5, 1e6, 1e7),
    labels = scales::comma_format()
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1200),
    breaks = seq(0, 1250, 250)
  ) +
  theme_classic(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  xlab("years in past") +
  ylab(expression(paste("N"[e], "*10"^3))) +
  annotation_logticks(sides = "b")


#####ROH ANALYSIS#####

cryege10000 <- read.csv("~/pulcher/het/cryege_v2024.1.fa.hifireads.10000_pi.txt", sep="\t")
cryege10000 <- subset(cryege10000, no_sites >= 5000)
crypul10000 <- read.csv("~/pulcher/het/crypul_v2024.1.fa.hifireads.10000_pi.txt", sep="\t")
crypul10000 <- subset(crypul10000, no_sites >= 5000)
crypul_omni10000 <- read.csv("~/pulcher/het/crypul_v2024.1.fa.omniCreads.10000_pi.txt", sep="\t")
crypul_omni10000 <- subset(crypul_omni10000, no_sites >= 5000)


smooth_ROH <- function(input_df) {
  count_diffs_smooth <- c()
  for (i in 1:nrow(input_df)) {
    if (i == 1) {
      count_diffs_smooth <- c(count_diffs_smooth, input_df$count_diffs[i])
    } else if (i == nrow(input_df)) {
      count_diffs_smooth <- c(count_diffs_smooth, input_df$count_diffs[i])
    } else {
      if (!is.na(input_df$count_diffs[i]) && !is.na(input_df$count_diffs[i-1]) && !is.na(input_df$count_diffs[i+1])) {
        if (input_df$count_diffs[i] > 2 && input_df$count_diffs[i-1] <= 2 && input_df$count_diffs[i+1] <= 2) {
          count_diffs_smooth <- c(count_diffs_smooth, 0)
        } else {
          count_diffs_smooth <- c(count_diffs_smooth, input_df$count_diffs[i])
        }
      } else {
        count_diffs_smooth <- c(count_diffs_smooth, input_df$count_diffs[i])
      }
    }
  }
  output_df <- cbind(input_df, count_diffs_smooth = as.data.frame(count_diffs_smooth))
  return(output_df)
}

find_ROH <- function(input_df) {
  output_df <- NULL
  begin <- c()
  end <- c()
  chroms <- c()
  nbins <- c()
  
  for (current_chrom in unique(input_df$chromosome)) {
    subset_input_df <- subset(input_df, chromosome == as.character(current_chrom))
    
    for (i in 1:nrow(subset_input_df)) {
      if (nrow(subset_input_df) == 1) {
        if (subset_input_df$count_diffs_smooth[i] <= 2) {
          start <- i
          chroms <- c(chroms, subset_input_df$chromosome[i])
          begin <- c(begin, subset_input_df$window_pos_1[i])
          end <- c(end, subset_input_df$window_pos_2[i])
          finish <- i + 1
          count <- finish - start
          nbins <- c(nbins, count)
        }
      } else {
        if (i == 1) {
          if (subset_input_df$count_diffs_smooth[i] <= 2) {
            chroms <- c(chroms, subset_input_df$chromosome[i])
            begin <- c(begin, subset_input_df$window_pos_1[i])
            start <- i
            if (subset_input_df$count_diffs_smooth[i + 1] > 2) {
              end <- c(end, subset_input_df$window_pos_2[i])
              finish <- i + 1
              count <- finish - start
              nbins <- c(nbins, count)
            }
          }
        }
        if (i > 1) {
          if (subset_input_df$count_diffs_smooth[i - 1] > 2 & subset_input_df$count_diffs_smooth[i] <= 2) {
            chroms <- c(chroms, subset_input_df$chromosome[i])
            begin <- c(begin, subset_input_df$window_pos_1[i])
            start <- i
          }
          if (i < nrow(subset_input_df)) {
            if (subset_input_df$count_diffs_smooth[i + 1] > 2 & subset_input_df$count_diffs_smooth[i] <= 2) {
              end <- c(end, subset_input_df$window_pos_2[i])
              finish <- i + 1
              count <- finish - start
              nbins <- c(nbins, count)
            }
          }
          if (i == nrow(subset_input_df)) {
            if (subset_input_df$count_diffs_smooth[nrow(subset_input_df)] <= 2) {
              end <- c(end, subset_input_df$window_pos_2[nrow(subset_input_df)])
              finish <- i + 1
              count <- finish - start
              nbins <- c(nbins, count)
            }
          }
        }
      }
    }
  }
  
  if (length(chroms) == length(begin) && length(begin) == length(end) && length(end) == length(nbins)) {
    output_df <- data.frame(chrom = chroms, begin = begin, end = end, nbins = nbins)
    output_df$length_bp <- (output_df$end - output_df$begin)
    output_df$length_mb <- (output_df$end - output_df$begin) / 1e6
  } else {
    warning("Vectors have different lengths; please check the logic.")
  }
  
  return(output_df)
}

cryege10000_smooth <- smooth_ROH(cryege10000)
cryege10000_ROH <- find_ROH(cryege10000_smooth)
cryege10000_ROH <- subset(cryege10000_ROH, nbins >= 100 & chrom != "chr7")
cryege10000_ROH$length_cat <- cut(cryege10000_ROH$length_mb, breaks = c(1, 5, 10, 15, 100))


crypul10000_smooth <- smooth_ROH(crypul10000)
crypul10000_ROH <- find_ROH(crypul10000_smooth)
crypul10000_ROH <- subset(crypul10000_ROH, nbins >= 100 & chrom != "chr7")

crypul_omni10000_smooth <- smooth_ROH(crypul_omni10000)
crypul_omni10000_ROH <- find_ROH(crypul_omni10000_smooth)
crypul_omni10000_ROH <- subset(crypul_omni10000_ROH, nbins >= 100 & chrom != "chr7")

cryege10000_ROH$length_cat <- cut(cryege10000_ROH$length_mb, breaks = c(1, 5, 10, 15, 100))
crypul10000_ROH$length_cat <- cut(crypul10000_ROH$length_mb, breaks = c(1, 5, 10, 15, 100))
crypul_omni10000_ROH$length_cat <- cut(crypul_omni10000_ROH$length_mb, breaks = c(1, 5, 10, 15, 100))

cryege10000_ROH <- cryege10000_ROH %>% arrange(length_bp) %>% mutate(length_sum = cumsum(length_bp)) %>% mutate(species="C. egeriae (HiFi)") %>% arrange(length_bp)
crypul10000_ROH <- crypul10000_ROH %>% arrange(length_bp) %>% mutate(length_sum = cumsum(length_bp)) %>% mutate(species="C. pulcher (HiFi)") %>% arrange(length_bp)
crypul_omni10000_ROH <- crypul_omni10000_ROH %>% arrange(length_bp) %>% mutate(length_sum = cumsum(length_bp)) %>% mutate(species="C. pulcher (omni-C)") %>% arrange(length_bp)

skink_ROH <- rbind(crypul10000_ROH,crypul_omni10000_ROH,cryege10000_ROH)
levels(skink_ROH$length_cat) <- list("1-5 Mb"="(1,5]","5-10 Mb"="(5,10]","10-15 Mb"="(10,15]",">15 Mb"="(15,100]")

skink_ROH$length_cat <- factor(skink_ROH$length_cat, levels=c(">15 Mb","10-15 Mb","5-10 Mb","1-5 Mb"))
skink_ROH$species <- factor(skink_ROH$species, levels=c("C. pulcher (HiFi)","C. pulcher (omni-C)","C. egeriae (HiFi)"))

ROH_plot <- ggplot() +
  geom_col(data=skink_ROH, aes(x=species, y=length_mb, fill=length_cat), color="black", position="stack", lwd=0.25) +
  geom_col(data=subset(skink_ROH, species == "C. pulcher (HiFi)"), aes(x=1, y=length_mb), fill="lightgrey", color="black", position="stack", lwd=0.25) +
  geom_col(data=subset(skink_ROH, species == "C. pulcher (omni-C)"), aes(x=2, y=length_mb), fill="lightgrey", color="black", position="stack", lwd=0.25) +
  theme_classic(base_size = 8) +
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.3, 0.75),
        legend.direction = "horizontal") +
  xlab("") +
  ylab("cumulative ROH length (Mb)") +
  scale_fill_manual("", values=c("#043265", "#1D4776", "#4F7094", "#9AADC2"))



#popgen_figure
### --- Combine full figure ---
het_plot <- plot_grid(fig_hifi, fig_omni, fig_egeriae, nrow = 3, rel_heights = c(1, 1, 1))
ROH_PSMC <- plot_grid(ROH_plot, PSMC_plot, ncol=2, rel_widths= c(1,1.5), labels = c("(g)","(h)"), label_size = 10)
popgen_plot <- plot_grid(het_plot, ROH_PSMC, nrow = 2, rel_heights= c(1.25,1))
popgen_plot
ggsave("~/pulcher/figs/fig3_popgen_plot.pdf",
       popgen_plot, height = 169, width = 169, units = "mm")


###################
#### Figure S2 ####
###################


read_and_normalize_coverage <- function(infile) {
  df <- read.csv(infile, sep = "\t", header = FALSE)
  names(df) <- c("chrom", "start", "end", "coverage")
  df$bin_width <- df$end - df$start
  mean_cov <- weighted.mean(df$coverage, w = df$bin_width)
  df$rel_cov <- df$coverage / mean_cov
  df
}


plot_coverage_allchr_dual <- function(
    infile1, infile2, ylims, xlab_text, plot_title,
    color_vector = c("#272e34", "#A69D97"), alpha = 0.8,
    label1 = "with Y", label2 = "no Y"
) {
  df1 <- read_and_normalize_coverage(infile1)
  df2 <- read_and_normalize_coverage(infile2)
  
  process_df <- function(df, label) {
    df %>%
      filter(bin_width == 1e6) %>%
      mutate(name = gsub("chr", "", chrom),
             name = factor(name, levels = c(as.character(1:15), "Y")),
             line_type = label)
  }
  
  df1_proc <- process_df(df1, label1)
  df2_proc <- process_df(df2, label2)
  df_combined <- bind_rows(df1_proc, df2_proc)
  
  contig_offsets <- df_combined %>%
    group_by(name) %>%
    summarise(chr_end = max(end), .groups = "drop") %>%
    mutate(chr_len_gap = chr_end + 1e7,
           offset = cumsum(chr_len_gap) - chr_len_gap) %>%
    select(name, offset)
  
  df_combined <- df_combined %>%
    left_join(contig_offsets, by = "name") %>%
    arrange(name, end) %>%
    mutate(pos2 = end + offset,
           pos = pos2 - bin_width)
  
  axis_df <- df_combined %>%
    group_by(name) %>%
    summarise(center = (max(pos2) + min(pos2)) / 2, .groups = "drop") %>%
    arrange(-center) %>%
    mutate(name_plot = str_replace_all(as.character(name), c('11'=' ', '13'=' ', '15'=' ')))
  
  ggplot(df_combined, aes(x = pos2 - 5e5, y = log2(rel_cov))) +
    geom_line(aes(color = chrom, linetype = line_type), alpha = alpha, size = 0.6) +
    scale_color_manual(values = rep(color_vector, length.out = length(unique(df_combined$chrom))), guide = "none") +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_x_continuous(labels = axis_df$name_plot, breaks = axis_df$center) +
    scale_y_continuous(limits = ylims, expand = c(0, 0)) +
    theme_classic(base_size = 8) +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top"
          ) +
    xlab(xlab_text) +
    ylab(expression(paste("log"[2], "(read-depth)", sep = ""))) +
    ggtitle(plot_title)
}

plot_coverage_singlechr_dual <- function(
    infile1, infile2, chromosome, ylims, xlab_text, plot_title,
    color = "#272e34", alpha = 1,
    label1 = "with Y", label2 = "no Y"
) {
  df1 <- read_and_normalize_coverage(infile1)
  df2 <- read_and_normalize_coverage(infile2)
  
  df1_chr <- df1 %>%
    filter(bin_width == 1e6, chrom == chromosome) %>%
    mutate(line_type = label1)
  
  df2_chr <- df2 %>%
    filter(bin_width == 1e6, chrom == chromosome) %>%
    mutate(line_type = label2)
  
  df_combined <- bind_rows(df1_chr, df2_chr)
  
  ggplot(df_combined, aes(x = start + 5e5, y = log2(rel_cov), linetype = line_type)) +
    geom_line(color = color, alpha = alpha, size = 0.6) +
    scale_x_continuous(breaks = seq(0, 70e6, 10e6), labels = seq(0, 70, 10)) +
    scale_y_continuous(limits = ylims, expand = c(0, 0)) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    theme_classic(base_size = 8) +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = "top"
          ) +
    xlab(xlab_text) +
    ylab(expression(paste("log"[2], "(read-depth)", sep = ""))) +
    ggtitle(plot_title)
}


crypul_WG <- plot_coverage_allchr_dual(
  infile1 = "~/pulcher/coverage/crypul_v2024.1.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  infile2 = "~/pulcher/coverage/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  ylims = c(-1.5, 1.5),
  xlab_text = "C. pulcher chromosomes",
  plot_title = "C. pulcher whole-genome coverage",
  color_vector = c("#272e34", "#A69D97"),
  alpha = 0.6,
  label1 = "scaf. 16 (putative Y) included",
  label2 = "scaf. 16 (putative Y) excluded"
)


crypul_chr7 <- plot_coverage_singlechr_dual(
  infile1 = "~/pulcher/coverage/crypul_v2024.1.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  infile2 = "~/pulcher/coverage/crypul_v2024.1_rmY.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  chromosome = "chr7",
  ylims = c(-1.5, 1.5),
  xlab_text = "C. pulcher chromosome 7 position (Mb)",
  plot_title = "C. pulcher coverage on chr7",
  color = "#272e34",
  alpha = 0.6,
  label1 = "scaf. 16 (putative Y) included",
  label2 = "scaf. 16 (putative Y) excluded"
)


cryege_WG <- plot_coverage_allchr_dual(
  infile1 = "~/pulcher/coverage/cryege_v2024.1.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  infile2 = "~/pulcher/coverage/cryege_v2024.1_rmY.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  ylims = c(-1.5, 1.5),
  xlab_text = "C. egeriae chromosomes",
  plot_title = "C. egeriae whole-genome coverage",
  color_vector = c("#092d44","#69b9cd"),
  alpha = 0.6,
  label1 = "scaf. 16 (putative Y) included",
  label2 = "scaf. 16 (putative Y) excluded"
)


cryege_chr7 <- plot_coverage_singlechr_dual(
  infile1 = "~/pulcher/coverage/cryege_v2024.1.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  infile2 = "~/pulcher/coverage/cryege_v2024.1_rmY.fa.hifireads.sorted.q20.bam_whole-genome-coverage_windows_1000000.bed",
  chromosome = "chr7",
  ylims = c(-1.5, 1.5),
  xlab_text = "C. egeriae chromosome 7 position (Mb)",
  plot_title = "C. egeriae coverage on chr7",
  color = "#092d44",
  alpha = 0.6,
  label1 = "scaf. 16 (putatuve Y) included",
  label2 = "scaf. 16 (putative Y) excluded"
)


coverage_WG <- plot_grid(crypul_WG, cryege_WG, ncol=1, labels = c("(a)","(c)"), label_size = 10)
coverage_chr7 <- plot_grid(crypul_chr7, cryege_chr7, ncol=1, labels = c("(b)","(d)"), label_size = 10)
coverage_fig <- plot_grid(coverage_WG, coverage_chr7, ncol=2, rel_widths = c(3,3))

ggsave("~/pulcher/figs/FigS2_coverage.pdf",
       coverage_fig, height = 100, width = 169, units = "mm")



###################
#### Figure S4 ####
###################

PSMC_hifi <- read.csv("~/pulcher/psmc/crypul_v2024.1.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.0.txt", header=FALSE, sep = "\t")
PSMC_omnic <- read.csv("~/pulcher/psmc/crypul_v2024.1.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.0.txt", header=FALSE, sep = "\t")
PSMC_hifi_callable <- read.csv("~/pulcher/psmc/crypul_v2024.1.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.0.txt", header=FALSE, sep = "\t")
PSMC_omnic_callable <- read.csv("~/pulcher/psmc/crypul_v2024.1.fa.omniCreads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.hifiomniC_callable.test.0.txt", header=FALSE, sep = "\t")
PSMC_egeriae <- read.csv("~/pulcher/psmc/cryege_v2024.1.fa.hifireads.sorted.q20.bam.autosomes.SNP_INVARIANT.FILT.PASS.0.txt", header=FALSE, sep = "\t")

PSMC_hifi$Source <- "C. pulcher (HiFi)"
PSMC_omnic$Source <- "C. pulcher (omni-C)"
PSMC_hifi_callable$Source <- "C. pulcher (HiFi - callable)"
PSMC_omnic_callable$Source <- "C. pulcher (omni-C - callable)"
PSMC_egeriae$Source <- "C. egeriae (HiFi)"

PSMC_combined <- rbind(PSMC_hifi, PSMC_omnic, PSMC_hifi_callable, PSMC_omnic_callable, PSMC_egeriae)

PSMC_supp <- ggplot() +
  geom_step(data=PSMC_combined, aes(x=V1, y=V2*10, color=Source)) +
  
  scale_x_log10(breaks = c(1e4,1e5,1e6, 1e7, 1e8),
                labels = c("1e4","1e5","1e6", "1e7", "1e8")) +
  scale_y_continuous(
    expand=c(0,0), limits = c(0,1500), breaks = c(seq(0,2000,250))) +
  theme_classic(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("years in past") +
  ylab(expression(paste("N"[e], "*", 10^{3},sep=" "))) +
  annotation_logticks(side = "b")

ggsave("~/pulcher/figs/FigS4_PSMC_supp.pdf",
       PSMC_supp, height = 115, width = 115, units = "mm")
