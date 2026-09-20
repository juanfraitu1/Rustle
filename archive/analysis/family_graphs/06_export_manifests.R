source("02_build_graphs.R")

DATA_DIR <- "data"

rbmy_vg <- readRDS(file.path(DATA_DIR, "rbmy_vg_seq.rds"))
export_family_manifest(rbmy_vg$exon_df, "RBMY",
                       file.path(DATA_DIR, "rbmy_manifest.tsv"))
message("Exported RBMY manifest: ", length(unique(rbmy_vg$exon_df$gene_id)),
        " loci → data/rbmy_manifest.tsv")

golga_vg <- readRDS(file.path(DATA_DIR, "golga_vg_seq.rds"))
export_family_manifest(golga_vg$exon_df, "GOLGA",
                       file.path(DATA_DIR, "golga_manifest.tsv"))
message("Exported GOLGA manifest: ", length(unique(golga_vg$exon_df$gene_id)),
        " loci → data/golga_manifest.tsv")

message("Done. Check data/*_manifest.tsv")
