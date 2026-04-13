# Build a custom OrgDb for Bombyx mori from NCBI GAF annotation
# Input:  data/bm_ncbi_gene_ontology.gaf
# Output: data/org.Bmori.eg.db_<version>.tar.gz (re-installed automatically)
library(AnnotationForge)
library(GO.db)
library(dplyr)

GAF_FILE  <- "data/bm_ncbi_gene_ontology.gaf"
TAXID     <- 7091L          # Bombyx mori NCBI taxonomy ID
GENUS     <- "Bombyx"
SPECIES   <- "mori"
VERSION   <- "0.3.0"    # R package version must be X.Y or X.Y.Z format
PKG_DIR   <- "data"

# ---- 1. Parse GAF ----
# Columns per GAF 2.2 spec (tab-separated, skip comment lines starting with !)
gaf_cols <- c("DB", "GeneID", "Symbol", "Qualifier", "GO_ID",
              "Reference", "Evidence", "With_From", "Aspect",
              "Gene_Name", "Synonym", "Type", "Taxon",
              "Date", "Assigned_By", "Annot_Ext", "Gene_Product_Form")

gaf <- read.table(GAF_FILE, sep = "\t", header = FALSE,
                  comment.char = "!", quote = "",
                  col.names = gaf_cols, fill = TRUE)

cat("GAF rows:", nrow(gaf), "\n")
cat("Unique genes:", length(unique(gaf$GeneID)), "\n")

# ---- 2. gene_info: GID + Symbol + Gene_Name ----
gene_info <- gaf |>
  select(GID = GeneID, SYMBOL = Symbol, GENENAME = Gene_Name) |>
  distinct(GID, .keep_all = TRUE) |>
  mutate(GID = as.character(GID))

cat("gene_info rows:", nrow(gene_info), "\n")

# ---- 3. gene2go: GID + GO + EVIDENCE ----
# makeOrgPackage requires exactly 3 columns: GID, GO, EVIDENCE
# Aspect filter: F=MF, P=BP, C=CC (exclude unrecognized codes)
gene2go <- gaf |>
  filter(GO_ID != "", Aspect %in% c("F", "P", "C")) |>
  mutate(
    GID      = as.character(GeneID),
    GO       = GO_ID,
    EVIDENCE = Evidence
  ) |>
  select(GID, GO, EVIDENCE) |>    # makeOrgPackage requires exactly these 3 cols
  distinct()

cat("gene2go rows:", nrow(gene2go), "\n")

# ---- 4. Make OrgDb via AnnotationForge ----
# makeOrgPackage expects a list of data frames; gene_info is mandatory
# Column GID must be "GID" as the primary key

pkg_name <- paste0("org.Bmori.eg.db")

# Remove existing package source directory before rebuild
unlink(file.path(PKG_DIR, pkg_name), recursive = TRUE)

makeOrgPackage(
  gene_info = gene_info,
  go        = gene2go,
  version   = VERSION,
  maintainer = "User <user@example.com>",
  author     = "User",
  outputDir  = PKG_DIR,
  tax_id     = as.character(TAXID),
  genus      = GENUS,
  species    = SPECIES,
  goTable    = "go",        # name of the go data frame passed above
  verbose    = TRUE
)

# ---- 5. Build and install tarball ----
pkg_path <- file.path(PKG_DIR, pkg_name)

# Remove previously built tarball to avoid confusion
old_tarballs <- list.files(PKG_DIR, pattern = paste0(pkg_name, "_.*\\.tar\\.gz"),
                           full.names = TRUE)
file.remove(old_tarballs)

# Build source tarball (outputs to current working directory)
build_out <- system2("R", args = c("CMD", "build", "--no-build-vignettes", pkg_path),
                     stdout = TRUE, stderr = TRUE)
cat(build_out, sep = "\n")

# Move tarball from cwd into data/
built <- list.files(".", pattern = paste0(pkg_name, "_.*\\.tar\\.gz"),
                    full.names = TRUE)
if (length(built) > 0) {
  tarball <- file.path(PKG_DIR, basename(built[1]))
  file.rename(built[1], tarball)
  cat("Tarball:", tarball, "\n")
  install.packages(tarball, repos = NULL, type = "source")
} else {
  # Fallback: install directly from source directory
  cat("Tarball build failed -- installing from source directory\n")
  install.packages(pkg_path, repos = NULL, type = "source")
}
cat("Installed:", pkg_name, "\n")

# ---- 6. Quick sanity check ----
# Reload to pick up new version
if (isNamespaceLoaded("org.Bmori.eg.db")) unloadNamespace("org.Bmori.eg.db")
library(org.Bmori.eg.db)
kt <- keytypes(org.Bmori.eg.db)
cat("Keytypes:", paste(kt, collapse = ", "), "\n")

# Use actual column names from keytypes
sym_col <- intersect(kt, c("SYMBOL", "Symbol"))[1]
test_gids <- head(keys(org.Bmori.eg.db, keytype = "GID"), 3)
test_res  <- AnnotationDbi::select(
  org.Bmori.eg.db,
  keys    = test_gids,
  keytype = "GID",
  columns = c(sym_col, "GO", "ONTOLOGY")
)
cat("Sample annotation:\n")
print(head(test_res, 6))

# Check coverage for our sig genes
sig_df  <- read.table("data/output/sig_splicing_genes_annotated.tsv",
                      header = TRUE, sep = "\t")
matched <- sum(as.character(sig_df$gene_id) %in% keys(org.Bmori.eg.db, keytype = "GID"))
cat("\nSig gene coverage:", matched, "/", nrow(sig_df), "\n")
