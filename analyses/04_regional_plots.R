# 4. regional plots
# use SNCA as example. 

#!/usr/bin/env Rscript

## ================================================================
## LocusZoom-style regional plot with gene track (SNCA example)
## - Uses PLINK 1.9 LD output and PLINK2 summary statistics
## - Build: hg38
## - Output: PNG/PDF + merged regional table
## ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(grid)
})

## ---------------- USER CONFIGURATION ----------------------------

# GWAS summary statistics (PLINK2 format)
sumstats_file <- "imputed_logistic_PLINK2_frequency.csv"
sumstats_sep  <- ","          # use "\t" for TSV

# PLINK LD output file (produced e.g. by commands below)
# Example:
#   plink2 --bfile /path/to/chr4_TWB_imputed \
#          --keep control_list.txt \
#          --make-bed --out chr4_control
#
#   plink19 --bfile chr4_control \
#           --ld-snp chr4:89759478:A:ATGCATATT \
#           --ld-window-kb 500 \
#           --ld-window 999999 \
#           --ld-window-r2 0 \
#           --r2 --out SNCA
ld_file    <- "SNCA.ld"
ld_prefix  <- "SNCA"          # used only in comments above

# GTF annotation for gene track (hg38)
gtf_path   <- "/Users/chuyungtsai/Desktop/Bioinfo/hg38/Homo_sapiens.GRCh38.109.gtf.gz"

# Region of interest: SNCA example
build      <- "hg38"
lead_chr   <- "4"
lead_bp    <- 89759478L
lead_var   <- "chr4:89759478:A:ATGCATATT"   # PLINK chr:BP:REF:ALT ID
window_kb  <- 500L                          # +/- 500 kb

# Output prefix
out_prefix <- "20251010_chr4_SNCA_region"

## ---------------- LOAD LD TABLE --------------------------------

ld <- fread(ld_file, data.table = FALSE)

need_ld_cols <- c("SNP_A", "SNP_B", "R2")
miss_ld <- setdiff(need_ld_cols, names(ld))
if (length(miss_ld)) {
  stop("LD table missing columns: ", paste(miss_ld, collapse = ", "))
}

# Keep only LD vs. the designated lead variant
ld <- ld %>%
  filter(SNP_A == lead_var) %>%
  transmute(VARID = SNP_B, R2 = as.numeric(R2))

## ---------------- LOAD SUMMARY STATISTICS ----------------------

ss <- fread(sumstats_file, sep = sumstats_sep, data.table = FALSE)

need <- c("CHR", "BP", "ID", "REF", "ALT", "A1", "P")
miss <- setdiff(need, names(ss))
if (length(miss)) {
  stop("Summary stats missing columns: ", paste(miss, collapse = ", "))
}

# Harmonize IDs:
# - CHR_chr: "chr4" style
# - uppercase alleles
# - VARID assumed to match PLINK ID (e.g. chr4:89759478:A:ATGCATATT)
ss <- ss %>%
  mutate(
    CHR_chr = ifelse(startsWith(as.character(CHR), "chr"),
                     as.character(CHR),
                     paste0("chr", CHR)),
    REF   = toupper(REF),
    ALT   = toupper(ALT),
    A1    = toupper(A1),
    VARID = ID
  )

## ---------------- REGION FILTER (± window_kb) -------------------

start_bp <- max(1L, lead_bp - window_kb * 1000L)
end_bp   <- lead_bp + window_kb * 1000L

ss_reg <- ss %>%
  filter(
    gsub("^chr", "", CHR_chr) == lead_chr,
    BP >= start_bp,
    BP <= end_bp
  )

## ---------------- MERGE LD ONTO SUMSTATS -----------------------

dat <- ss_reg %>%
  left_join(ld, by = "VARID")

# Fallback: if most R2 are NA, try flipped allele order(chr:BP:ALT:REF)
need_fallback <- sum(!is.na(dat$P) & is.na(dat$R2)) >
  0.5 * sum(!is.na(dat$P))

if (need_fallback) {
  message("Many LD matches are missing; trying fallback join with flipped alleles...")
  ss_reg_fallback <- ss_reg %>%
    mutate(VARID_FLIP = paste(CHR_chr, BP, ALT, REF, sep = ":"))
  ld_fallback <- ld %>%
    rename(VARID_FLIP = VARID)
  
  dat2 <- ss_reg_fallback %>%
    left_join(ld_fallback, by = "VARID_FLIP")
  
  dat$R2 <- ifelse(is.na(dat$R2) & !is.na(dat2$R2), dat2$R2, dat$R2)
}

## ---------------- ENSURE LEAD VARIANT IS PRESENT ---------------

lead_row <- ss %>% filter(VARID == lead_var)
if (nrow(lead_row) == 0L) {
  # placeholder if lead variant is not in summary stats
  lead_row <- data.frame(
    CHR_chr = paste0("chr", lead_chr),
    BP      = lead_bp,
    ID      = NA_character_,
    REF     = "A",
    ALT     = "ATGCATATT",
    A1      = "ATGCATATT",
    P       = NA_real_,
    VARID   = lead_var,
    stringsAsFactors = FALSE
  )
}
lead_row$R2 <- 1.0

dat <- dat %>%
  bind_rows(lead_row %>% anti_join(dat, by = "VARID")) %>%
  mutate(
    neglog10P = ifelse(!is.na(P) & P > 0, -log10(P), NA_real_),
    is_lead   = VARID == lead_var
  ) %>%
  mutate(
    ld_bin = cut(
      R2,
      breaks = c(-0.001, 0.2, 0.4, 0.6, 0.8, 0.999999, 1.01),
      labels = c("<0.2", "0.2–0.4", "0.4–0.6",
                 "0.6–0.8", "0.8–1.0", "Lead"),
      include.lowest = TRUE
    ),
    ld_bin = ifelse(is_lead, "Lead", as.character(ld_bin)),
    ld_bin = factor(ld_bin,
                    levels = c("<0.2", "0.2–0.4", "0.4–0.6",
                               "0.6–0.8", "0.8–1.0", "Lead"))
  )

## ---------------- LABELS (rs5860181 + top hit) -----------------

# Force label mapping for specific variants
custom_labels <- c(
  "chr4:89759478:A:ATGCATATT" = "rs5860181"
)

dat <- dat %>%
  mutate(label_override = custom_labels[VARID])

lab_top <- dat %>%
  filter(!is.na(P)) %>%
  arrange(P) %>%
  slice_head(n = 1L)

lab_lead <- dat %>%
  filter(is_lead)

lab <- bind_rows(lab_top, lab_lead) %>%
  distinct(VARID, .keep_all = TRUE) %>%
  mutate(
    lbl = coalesce(
      label_override,
      if_else(!is.na(ID) & ID != ".", ID, NA_character_),
      VARID
    )
  )

## ---------------- GENE TRACK (GTF) -----------------------------

p_gene <- NULL  # default if GTF/genes are not available

if (file.exists(gtf_path)) {
  gtf <- fread(
    gtf_path,
    sep         = "\t",
    header      = FALSE,
    quote       = "",
    na.strings  = c("", "NA")
  )
  
  setnames(
    gtf,
    paste0("V", 1:9),
    c("seqname", "source", "feature", "start", "end",
      "score", "strand", "frame", "attr")
  )
  
  # Limit to region and chromosome
  gtf <- gtf %>%
    filter(
      gsub("^chr", "", seqname) == as.character(lead_chr),
      start <= end_bp,
      end   >= start_bp
    )
  
  genes_info <- gtf %>%
    filter(feature == "gene") %>%
    transmute(
      gene_id   = str_match(attr, 'gene_id "([^"]+)"')[, 2],
      gene_name = str_match(attr, 'gene_name "([^"]+)"')[, 2],
      gene_start = pmax(start, start_bp),
      gene_end   = pmin(end,   end_bp),
      strand
    ) %>%
    filter(!is.na(gene_id))
  
  exons_tbl <- gtf %>%
    filter(feature == "exon") %>%
    transmute(
      gene_id       = str_match(attr, 'gene_id "([^"]+)"')[, 2],
      transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[, 2],
      start, end, strand
    ) %>%
    filter(!is.na(gene_id), !is.na(transcript_id))
  
  if (nrow(genes_info) > 0L && nrow(exons_tbl) > 0L) {
    # Choose canonical transcript (longest exonic span in window)
    tx_len <- exons_tbl %>%
      mutate(exon_len = pmin(end, end_bp) - pmax(start, start_bp) + 1L) %>%
      group_by(gene_id, transcript_id) %>%
      summarise(tx_exonic_len = sum(pmax(exon_len, 0L)), .groups = "drop")
    
    canon_tx <- tx_len %>%
      group_by(gene_id) %>%
      slice_max(tx_exonic_len, n = 1L, with_ties = FALSE) %>%
      ungroup()
    
    # Canonical exons with gene_name
    exons_canon <- exons_tbl %>%
      semi_join(canon_tx, by = c("gene_id", "transcript_id")) %>%
      left_join(
        genes_info %>% select(gene_id, gene_name, strand_gene = strand),
        by = "gene_id"
      ) %>%
      mutate(
        strand = if_else(!is.na(strand), strand, strand_gene),
        xmin   = pmax(start, start_bp),
        xmax   = pmin(end,   end_bp)
      ) %>%
      filter(xmax > xmin) %>%
      select(gene_id, gene_name, transcript_id, strand, xmin, xmax)
    
    # Canonical introns
    introns_canon <- exons_canon %>%
      arrange(gene_id, transcript_id, xmin) %>%
      group_by(gene_id, transcript_id) %>%
      mutate(prev_xmax = lag(xmax)) %>%
      ungroup() %>%
      filter(!is.na(prev_xmax)) %>%
      transmute(
        gene_id, transcript_id, gene_name, strand,
        intron_start = prev_xmax,
        intron_end   = xmin
      ) %>%
      filter(intron_end > intron_start)
    
    # Pack genes into non-overlapping rows
    genes_packed <- genes_info %>%
      arrange(gene_start) %>%
      mutate(row = {
        rows_end <- numeric(0)
        out      <- integer(n())
        for (i in seq_len(n())) {
          placed <- FALSE
          for (r in seq_along(rows_end)) {
            if (gene_start[i] > rows_end[r] + 1e4) {
              out[i]     <- r
              rows_end[r] <- gene_end[i]
              placed <- TRUE
              break
            }
          }
          if (!placed) {
            rows_end   <- c(rows_end, gene_end[i])
            out[i]     <- length(rows_end)
          }
        }
        out
      })
    
    # Attach rows
    exons_canon <- exons_canon %>%
      left_join(genes_packed %>% select(gene_id, row), by = "gene_id") %>%
      mutate(y = max(genes_packed$row, na.rm = TRUE) - row + 1)
    
    introns_canon <- introns_canon %>%
      left_join(genes_packed %>% select(gene_id, row), by = "gene_id") %>%
      mutate(
        y = max(genes_packed$row, na.rm = TRUE) - row + 1,
        draw_x    = if_else(strand == "+", intron_start, intron_end),
        draw_xend = if_else(strand == "+", intron_end,   intron_start)
      )
    
    # Clip to window
    exons_plot <- exons_canon %>%
      mutate(
        x0 = pmax(xmin, start_bp),
        x1 = pmin(xmax, end_bp)
      ) %>%
      filter(!is.na(gene_name), !is.na(y), x1 > x0)
    
    introns_plot <- introns_canon %>%
      mutate(
        x0 = pmax(draw_x,    start_bp),
        x1 = pmin(draw_xend, end_bp)
      ) %>%
      filter(!is.na(gene_name), !is.na(y), x1 > x0)
    
    if (nrow(exons_plot) > 0L || nrow(introns_plot) > 0L) {
      y_top <- max(c(exons_plot$y, introns_plot$y), na.rm = TRUE)
      
      p_gene <- ggplot() +
        geom_segment(
          data = introns_plot,
          aes(x = x0, xend = x1, y = y, yend = y),
          linewidth = 0.6,
          colour    = "grey20",
          arrow     = arrow(
            type   = "closed",
            length = unit(0.14, "cm")
          ),
          inherit.aes = FALSE
        ) +
        geom_rect(
          data = exons_plot,
          aes(
            xmin = x0, xmax = x1,
            ymin = y - 0.18, ymax = y + 0.18
          ),
          fill        = "grey30",
          colour      = "black",
          linewidth   = 0.2,
          inherit.aes = FALSE
        ) +
        geom_text(
          data = genes_packed %>%
            mutate(y = max(row, na.rm = TRUE) - row + 1),
          aes(x = (gene_start + gene_end) / 2, y = y - 0.35,
              label = gene_name),
          size        = 3,
          vjust       = 1,
          lineheight  = 0.9,
          inherit.aes = FALSE
        ) +
        scale_x_continuous(
          limits = c(start_bp, end_bp),
          labels = label_number(big.mark = ",")
        ) +
        scale_y_continuous(
          NULL, breaks = NULL,
          limits = c(0.4, y_top + 0.8),
          expand = c(0, 0)
        ) +
        labs(x = paste0("Position on chr", lead_chr, " (", build, ")")) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "black"),
          axis.title.y     = element_blank(),
          plot.margin      = margin(0, 12, 8, 12)
        ) +
        coord_cartesian(clip = "off")
    }
  }
}

## ---------------- MAIN LocusZoom-STYLE PANEL -------------------

lz_cols <- c(
  "<0.2"    = "#313695",  # dark blue
  "0.2–0.4" = "#74ADD1",  # light blue
  "0.4–0.6" = "#66BD63",  # green
  "0.6–0.8" = "#FDAE61",  # orange
  "0.8–1.0" = "#D73027",  # red
  "Lead"    = "#5E3C99"   # purple for lead
)

p_main <- ggplot(dat, aes(x = BP, y = neglog10P)) +
  geom_point(
    aes(fill = ld_bin, shape = is_lead),
    size   = 2.6,
    alpha  = 0.95,
    stroke = 0.6,
    colour = "black"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 23, `FALSE` = 21),   # diamond vs circle
    guide  = "none"
  ) +
  scale_fill_manual(
    values    = lz_cols,
    drop      = FALSE,
    na.value  = "grey70",
    name      = expression(r^2)
  ) +
  scale_x_continuous(
    limits = c(start_bp, end_bp),
    labels = label_number(big.mark = ",")
  ) +
  scale_y_continuous(
    name   = expression(-log[10](italic(P))),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    x = paste0("Position on chr", lead_chr, " (", build, ")"),
    title    = NULL,
    subtitle = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = c(0.12, 0.78),
    legend.background = element_rect(color = "black", fill = "white"),
    legend.margin     = margin(4, 6, 4, 6),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    panel.border      = element_rect(colour = "black"),
    plot.margin       = margin(6, 12, 2, 12)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21, size = 3, colour = "black", alpha = 1),
      ncol         = 1,
      byrow        = TRUE
    )
  )

if (nrow(lab) > 0L) {
  p_main <- p_main +
    ggrepel::geom_text_repel(
      data    = lab,
      aes(x = BP, y = neglog10P, label = lbl),
      size    = 3.2,
      min.segment.length = 0,
      box.padding        = 0.35,
      point.padding      = 0.25,
      max.overlaps       = 100,
      segment.size       = 0.25
    )
}

## ---------------- COMBINE PANELS & SAVE ------------------------

if (!is.null(p_gene)) {
  p_final   <- p_main / p_gene + plot_layout(heights = c(4, 1))
  height_out <- 6.2
} else {
  p_final   <- p_main
  height_out <- 5.5
}

ggsave(
  paste0(out_prefix, "_locus.png"),
  p_final,
  width  = 8.0,
  height = height_out,
  dpi    = 350
)

ggsave(
  paste0(out_prefix, "_locus.pdf"),
  p_final,
  width  = 8.0,
  height = height_out
)
