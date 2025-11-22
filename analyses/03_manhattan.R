# 3. Manhattan Plot
# ===== Manhattan plot (publication-ready, JPG output) =====
# Dependencies
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid) 

# ---- I/O + basics ----
infile  <- "sumstats.csv"                      # your GWAS summary stats
outfile <- "Figure1_Manhattan_highlighted.jpg"

gws_thr  <- 5e-8
sugg_thr <- 1e-5

# ---- loci to highlight (fill/adjust positions as needed; GRCh38) ----
# ---- loci to highlight (GRCh38) ----
# WINDOW_KB controls half-width of the orange band around each locus.
targets <- data.table(
  GENE = c("SNCA","MCCC1","GCH1","PPARGC1A","GALNT13","VPS13C","SCARB2","LRRK2"),
  CHR  = c(4,       3,       14,    4,         2,         15,       4,        12),
  POS  = c(89759478,183028420,54885803,24009221,154088802,61994009,76217199,40363526),
  WINDOW_KB = c(300, 300, 300, 300, 300, 300, 300, 300)
)


# ---- load & clean ----
ss <- fread(infile)
names(ss) <- toupper(names(ss))
if ("POS" %in% names(ss)) ss[, BP := POS]
stopifnot(all(c("CHR","BP","P") %in% names(ss)))

# normalize chromosomes 1..22
ss[, CHR := gsub("^chr", "", CHR, ignore.case = TRUE)]
suppressWarnings(ss[, CHR := as.integer(CHR)])
ss <- ss[CHR %in% 1:22]
# valid positions and p-values
ss <- ss[is.finite(BP) & is.finite(P) & P > 0 & P <= 1]
ss[P == 0, P := .Machine$double.xmin]
ss[, LOGP := -log10(P)]
ss[LOGP > 300, LOGP := 300]   # optional cap

# ---- cumulative coords & axis ----
setorder(ss, CHR, BP)
chroms <- sort(unique(ss$CHR))
chr_len <- ss[, .(CHR_LEN = max(BP, na.rm = TRUE)), by = CHR][data.table(CHR = chroms), on="CHR"]
chr_len[, CHR_LEN := as.numeric(CHR_LEN)]
chr_len[, CHR_OFFSET := c(0, head(cumsum(as.numeric(CHR_LEN)), -1))]
ss <- merge(ss, chr_len, by = "CHR", all.x = TRUE, sort = FALSE)
ss[, BP_CUM := BP + CHR_OFFSET]
axis_df <- chr_len[, .(CHR = chroms, center = CHR_OFFSET + CHR_LEN/2)]

# ---- color scheme: alternating grey / light-blue background ----
base_cols <- c("#b3b3b3", "#89c9f9")
alt_cols  <- setNames(rep(base_cols, length.out = length(chroms)), chroms)

# ---- build highlight mask for all provided loci ----
# keep only rows with valid CHR & POS
t_ok <- targets[is.finite(CHR) & is.finite(POS)]
if (nrow(t_ok) < nrow(targets)) {
  skipped <- setdiff(targets$GENE, t_ok$GENE)
  if (length(skipped)) message("Skipping highlights with missing CHR/POS: ",
                               paste(skipped, collapse=", "))
}

# tag highlighted variants
ss[, HILITE := FALSE]
if (nrow(t_ok)) {
  for (i in seq_len(nrow(t_ok))) {
    r <- t_ok[i]
    lo <- r$POS - r$WINDOW_KB*1000
    hi <- r$POS + r$WINDOW_KB*1000
    ss[CHR == r$CHR & BP >= lo & BP <= hi, HILITE := TRUE]
  }
}

# ---- choose labels: ONE per highlighted region (lead SNP) ----
labels <- data.table()
if (nrow(t_ok)) {
  for (i in seq_len(nrow(t_ok))) {
    r  <- t_ok[i]
    lo <- r$POS - r$WINDOW_KB*1000
    hi <- r$POS + r$WINDOW_KB*1000
    d  <- ss[CHR == r$CHR & BP >= lo & BP <= hi]
    if (nrow(d)) {
      # pick the single lead SNP (smallest P, i.e., largest -log10P)
      lead <- d[which.min(P)]
      lead[, GENE := r$GENE]
      labels <- rbind(labels, lead, fill = TRUE)
    }
  }
}

# text to show: prefer SNP, else gene name, else chr:bp
if ("SNP" %in% names(labels)) {
  labels[, LABEL := fifelse(!is.na(SNP) & SNP != "", SNP, GENE)]
} else {
  labels[, LABEL := GENE]
}
labels[is.na(LABEL) | LABEL == "", LABEL := paste0("chr", CHR, ":", BP)]


# ---- plot ----
p <- ggplot() +
  # background
  geom_point(data = ss[HILITE == FALSE],
             aes(BP_CUM, LOGP, color = factor(CHR)),
             size = 1.1, alpha = 0.65, stroke = 0) +
  scale_color_manual(values = alt_cols, guide = "none") +
  # all highlights in orange
  geom_point(data = ss[HILITE == TRUE],
             aes(BP_CUM, LOGP),
             size = 1.6, alpha = 0.95, color = "#F28E2B") +
  # optional threshold lines (uncomment if desired)
  # geom_hline(yintercept = -log10(gws_thr), color = "firebrick2", linetype = "dashed") +
  # geom_hline(yintercept = -log10(sugg_thr), color = "steelblue3", linetype = "dotted") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  labs(x = "BPcum", y = expression(-log[10](P))) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6))
  )

# labels: rounded white boxes with black border + leader lines
if (nrow(labels)) {
  p <- p +
    geom_point(data = labels, aes(BP_CUM, LOGP), size = 1.6, color = "#F28E2B") +
    ggrepel::geom_label_repel(
      data = labels,
      aes(BP_CUM, LOGP, label = LABEL),
      size = 3.1, label.size = 0.4, label.r = unit(0.18, "lines"),
      label.padding = unit(0.16, "lines"),
      color = "black", fill = "white", segment.color = "black",
      min.segment.length = 0, box.padding = 0.25, point.padding = 0.22,
      max.overlaps = Inf, force = 1
    )
}

jpeg(outfile, width = 2400, height = 1200, res = 220)
print(p)
dev.off()

cat("Saved:", normalizePath(outfile), "\n")
