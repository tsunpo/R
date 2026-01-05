# ============================================================
# Classify PacBio Iso-Seq reads around exon 10: Include10 / Skip10 / IR9
# Version-agnostic (no junctions(..., with.revmap=TRUE) dependency)
# Requires: Rsamtools, GenomicAlignments, GenomicRanges, dplyr, tibble, readr
# ------------------------------------------------------------
# Locus (hg38 example for CBL; change to your coordinates)
gene         <- "CBL"
chrom        <- "chr11"
strand       <- "-"            # gene strand (for interpretation only)
exon9_end    <- 119284963      # last base of exon 9
exon10_start <- 119284969      # first base of exon 10
exon10_end   <- 119285186      # last base of exon 10
exon11_start <- 119285559      # first base of exon 11

bam_wt  <- "/lustre/scratch125/casm/staging/team294/ty2/SSC/analysis/ssc/plots/CBL/11_119284966_T_A/origUM/BAM/WT.splice.bam"
bam_mut <- "/lustre/scratch125/casm/staging/team294/ty2/SSC/analysis/ssc/plots/CBL/11_119284966_T_A/origUM/BAM/MUT.splice.bam"

min_mapq  <- 20
junc_fuzz <- 3   # tolerance (bp) when matching donor/acceptor breakpoints
flank     <- 200 # window around locus to read from BAM
# ============================================================

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicAlignments)
	library(GenomicRanges)
	library(S4Vectors)
	library(BiocGenerics)
	library(dplyr); library(tibble); library(readr); library(stringr)
})

# ---- Region of interest -------------------------------------------------------
left  <- min(exon9_end, exon10_start) - flank
right <- max(exon10_end, exon11_start) + flank
roi   <- GRanges(seqnames = chrom, ranges = IRanges(start = left, end = right))

param <- ScanBamParam(
	which = roi,
	what  = c("qname","flag","rname","strand","pos","qwidth","cigar","mapq")
)

read_bam_as_gr <- function(bam) {
	aln <- readGAlignments(bam, param = param, use.names = TRUE)
	aln <- aln[mcols(aln)$mapq >= min_mapq]
	aln
}

aln_wt  <- read_bam_as_gr(bam_wt)
aln_mut <- read_bam_as_gr(bam_mut)

# ---- Extract splice junctions from CIGAR (N ops), per-read -------------------
# This avoids junctions(..., with.revmap=TRUE) and works across Bioc versions.
# ---- Extract splice-like gaps from CIGAR: N + large D ------------------------
junctions_from_aln <- function(aln, D_MIN = 20L) {
	if (length(aln) == 0) {
		return(tibble::tibble(qname = character(), jstart = integer(), jend = integer()))
	}
	cig <- GenomicAlignments::cigar(aln)
	aln_starts <- BiocGenerics::start(aln)
	qn <- names(aln); if (is.null(qn)) qn <- as.character(seq_along(aln))
	
	# 1) True splicing gaps (N)
	n_ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(
		cig, ops = "N", reduce.ranges = FALSE, drop.empty.ranges = TRUE
	)
	n_per_read <- S4Vectors::elementNROWS(n_ranges)
	n_idx_rep  <- rep.int(seq_along(aln), n_per_read)
	n_flat     <- unlist(n_ranges, use.names = FALSE)
	
	jt_n <- if (length(n_flat) > 0) tibble::tibble(
		qname  = qn[n_idx_rep],
		jstart = BiocGenerics::start(n_flat) + aln_starts[n_idx_rep] - 1L,
		jend   = BiocGenerics::end(n_flat)   + aln_starts[n_idx_rep] - 1L
	) else tibble::tibble(qname=character(), jstart=integer(), jend=integer())
	
	# 2) Large deletions (D) used as introns by DNA presets
	d_ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(
		cig, ops = "D", reduce.ranges = FALSE, drop.empty.ranges = TRUE
	)
	d_per_read <- S4Vectors::elementNROWS(d_ranges)
	d_idx_rep  <- rep.int(seq_along(aln), d_per_read)
	d_flat     <- unlist(d_ranges, use.names = FALSE)
	
	jt_d <- if (length(d_flat) > 0) {
		keep <- BiocGenerics::width(d_flat) >= D_MIN
		if (any(keep)) {
			tibble::tibble(
				qname  = qn[d_idx_rep][keep],
				jstart = BiocGenerics::start(d_flat)[keep] + aln_starts[d_idx_rep][keep] - 1L,
				jend   = BiocGenerics::end(d_flat)[keep]   + aln_starts[d_idx_rep][keep] - 1L
			)
		} else tibble::tibble(qname=character(), jstart=integer(), jend=integer())
	} else tibble::tibble(qname=character(), jstart=integer(), jend=integer())
	
	dplyr::bind_rows(jt_n, jt_d)
}

j_wt_tbl  <- junctions_from_aln(aln_wt)
j_mut_tbl <- junctions_from_aln(aln_mut)

# ---- Helpers -----------------------------------------------------------------
near <- function(x, target, fuzz) abs(x - target) <= fuzz

# Junction breakpoints (genomic)
j_9_10_donor   <- exon9_end
j_9_10_accept  <- exon10_start
j_10_11_donor  <- exon10_end
j_10_11_accept <- exon11_start
j_9_11_donor   <- exon9_end
j_9_11_accept  <- exon11_start

# ---- Collapse junctions into boolean flags per read --------------------------
per_read_flags <- function(junc_tbl, aln, junc_fuzz = 3) {
	base <- tibble(
		qname = names(aln) %||% as.character(seq_along(aln)),
		has_9_10  = FALSE,
		has_10_11 = FALSE,
		has_9_11  = FALSE
	)
	
	if (nrow(junc_tbl) > 0) {
		jdf <- junc_tbl %>%
			mutate(
				is_9_10  = near(jstart, j_9_10_donor,   junc_fuzz) & near(jend, j_9_10_accept,  junc_fuzz),
				is_10_11 = near(jstart, j_10_11_donor,  junc_fuzz) & near(jend, j_10_11_accept, junc_fuzz),
				is_9_11  = near(jstart, j_9_11_donor,   junc_fuzz) & near(jend, j_9_11_accept,  junc_fuzz)
			) %>%
			group_by(qname) %>%
			summarise(
				has_9_10  = any(is_9_10),
				has_10_11 = any(is_10_11),
				has_9_11  = any(is_9_11),
				.groups = "drop"
			)
		
		base <- base %>% left_join(jdf, by = "qname") %>%
			mutate(
				has_9_10  = coalesce(.data$has_9_10.y,  .data$has_9_10.x),
				has_10_11 = coalesce(.data$has_10_11.y, .data$has_10_11.x),
				has_9_11  = coalesce(.data$has_9_11.y,  .data$has_9_11.x)
			) %>%
			select(qname, has_9_10, has_10_11, has_9_11)
	}
	
	# IR9: read has a continuous aligned block that covers exon10_start (± fuzz)
	# WITHOUT having a 9→10 splice junction.
	bl <- GenomicAlignments::grglist(aln)   # exonic blocks per read (split by 'N')
	if (length(bl) == 0) {
		base$IR9 <- FALSE
		return(base)
	}
	bl_flat     <- unlist(bl, use.names = FALSE)
	block_read  <- rep(base$qname, elementNROWS(bl))
	
	acc_window <- GRanges(
		seqnames = chrom,
		ranges   = IRanges(exon10_start - junc_fuzz, exon10_start + junc_fuzz)
	)
	ov <- findOverlaps(acc_window, bl_flat, ignore.strand = TRUE)
	reads_covering_acceptor <- unique(block_read[subjectHits(ov)])
	
	base$IR9 <- base$qname %in% reads_covering_acceptor & !base$has_9_10
	base
}

flags_wt  <- per_read_flags(j_wt_tbl,  aln_wt,  junc_fuzz = junc_fuzz)
flags_mut <- per_read_flags(j_mut_tbl, aln_mut, junc_fuzz = junc_fuzz)

# ---- Final per-read classification -------------------------------------------
classify_reads <- function(flags_tbl) {
	flags_tbl %>%
		mutate(class = case_when(
			has_9_11 ~ "Skip10",
			has_9_10 & has_10_11 ~ "Include10",
			IR9 ~ "IR9",
			TRUE ~ "Other"
		))
}

wt_cls  <- classify_reads(flags_wt)  %>% mutate(sample = "WT")
mut_cls <- classify_reads(flags_mut) %>% mutate(sample = "MUT")
all_cls <- bind_rows(wt_cls, mut_cls)

# ---- Summaries (namespaced: no conflicts) ------------------------------------
counts <- all_cls %>%
	dplyr::count(sample, class, name = "n_reads") %>%
	dplyr::group_by(sample) %>%
	dplyr::mutate(frac = n_reads / sum(n_reads)) %>%
	dplyr::ungroup()

psi <- all_cls %>%
	dplyr::mutate(I = class == "Include10",
					  E = class == "Skip10") %>%
	dplyr::group_by(sample) %>%
	dplyr::summarise(
		I = sum(I),
		E = sum(E),
		PSI_exon10 = ifelse(I + E > 0, I / (I + E), NA_real_),
		.groups = "drop"
	)

readr::write_csv(counts, "exon10_class_counts.csv")
readr::write_csv(psi,    "exon10_PSI.csv")

print(counts)
print(psi)
