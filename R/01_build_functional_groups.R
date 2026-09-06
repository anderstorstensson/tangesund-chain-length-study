# =============================================================================
# Functional group, chain-forming flag, trophy and display name per IFCB class
# (WoRMS-backed via SHARK4R::assign_phytoplankton_group()).
#
# Reads:  data/SHARK_PlanktonImaging_*/received_data/class_names.txt (all deliveries),
#         data/roi_level/roi_2016.txt.gz; queries WoRMS unless data/.worms_assign_cache.rds
#         matches the query set.
# Writes: data/functional_group_map.csv
# Deps:   SHARK4R, data.table
# Run:    Rscript R/01_build_functional_groups.R
# =============================================================================
suppressPackageStartupMessages({ library(SHARK4R); library(data.table) })

cfg <- list(
  cn_glob  = "data/SHARK_PlanktonImaging_*/received_data/class_names.txt",
  roi_file = "data/roi_level/roi_2016.txt.gz",
  out_csv  = "data/functional_group_map.csv",
  chain_frac = 0.5      # class is 'chain-forming' if >=50% of its ROIs are YOLO chain-counted
)

# ---- union taxonomy across deliveries ---------------------------------------
cn <- unique(rbindlist(lapply(Sys.glob(cfg$cn_glob), function(f)
  fread(f, sep = "\t", colClasses = "character")), fill = TRUE), by = "class")
cn <- cn[class != "" & !is.na(class)]
cat(sprintf("union taxonomy: %d classes\n", nrow(cn)))

# ---- base group via SHARK4R (WoRMS), with the WoRMS class returned ------------
# query name: scientificname, else class_clean; unique names only, junk labels skipped
junklab <- c("Beads","Debris","Air_bubbles","Bubbles","mixed","Unicells","unclassified","Zooplankton")
cn[, qname := fifelse(!is.na(scientificname) & scientificname != "", scientificname, class_clean)]
qy <- unique(cn[!class %in% junklab & !is.na(qname) & qname != "",
                .(qname, aphia = suppressWarnings(as.integer(aphia_id)))], by = "qname")
# cache the WoRMS assignment (keyed on the query set) so re-runs are instant offline
wcache <- "data/.worms_assign_cache.rds"
key <- list(qname = sort(qy$qname))
if (file.exists(wcache) && identical(readRDS(wcache)$key, key)) {
  cat("WoRMS assignment cache hit\n"); g <- readRDS(wcache)$g
} else {
  cat(sprintf("querying WoRMS for %d unique taxon names (of %d classes)\n", nrow(qy), nrow(cn)))
  g <- assign_phytoplankton_group(scientific_names = qy$qname, aphia_ids = qy$aphia,
                                  return_class = TRUE, verbose = FALSE)
  setDT(g); saveRDS(list(key = key, g = g), wcache)
}
setDT(g)
cn[, worms_class_api := g$class[match(qname, g$scientific_name)]]
cn[, base_group     := g$plankton_group[match(qname, g$scientific_name)]]
# WoRMS class: API result, else the delivery's own worms_class
cn[, wclass := fifelse(!is.na(worms_class_api) & worms_class_api != "", worms_class_api, worms_class)]

# ---- refine the 'Other' bucket by WoRMS class, then by phylum fallback --------
ciliate  <- c("Litostomatea","Oligotrichea","Spirotrichea","Prostomatea","Heterotrichea","Karyorelictea","Nassophorea","Phyllopharyngea")
crypto   <- c("Cryptophyceae","Cryptophyta incertae sedis")
chloro   <- c("Pyramimonadophyceae","Chlorophyceae","Ulvophyceae","Trebouxiophyceae","Chlorodendrophyceae","Nephroselmidophyceae","Mamiellophyceae")
hapto    <- c("Coccolithophyceae","Prymnesiophyceae","Pavlovophyceae")
metazoan <- c("Bivalvia","Gastropoda","Copepoda","Maxillopoda","Appendicularia")
fine_of <- function(wc) {
  fifelse(wc %in% ciliate, "ciliate",
  fifelse(wc %in% crypto, "cryptophyte",
  fifelse(wc %in% chloro, "chlorophyte",
  fifelse(wc %in% hapto, "haptophyte",
  fifelse(wc == "Chrysophyceae", "chrysophyte",
  fifelse(wc == "Euglenophyceae", "euglenophyte",
  fifelse(wc == "Dictyochophyceae", "dictyochophyte",
  fifelse(wc %in% c("Raphidophyceae","Xanthophyceae"), "other stramenopile",
  fifelse(wc %in% c("Thecofilosea","Imbricatea","Katablepharidophyceae"), "heterotrophic flagellate",
  fifelse(wc %in% metazoan, "metazoan",
  "unassigned"))))))))))
}
cn[, group_fine := fcase(
  base_group == "Diatoms",         "diatom",
  base_group == "Dinoflagellates", "dinoflagellate",
  base_group == "Cyanobacteria",   "cyanobacteria",
  default = fine_of(wclass))]
# phylum fallback for taxa with no WoRMS class
cn[group_fine == "unassigned", group_fine := fcase(
  worms_phylum == "Ciliophora",                              "ciliate",
  worms_phylum == "Cryptophyta",                             "cryptophyte",
  worms_phylum %in% c("Haptophyta","Prymnesiophyta"),        "haptophyte",
  worms_phylum == "Cyanobacteria",                           "cyanobacteria",
  worms_phylum == "Chlorophyta",                             "chlorophyte",
  worms_phylum == "Euglenozoa",                              "euglenophyte",
  worms_phylum %in% c("Echinodermata","Mollusca","Arthropoda","Chordata","Cnidaria"), "metazoan",
  default = "unassigned")]
# non-plankton and junk classifier labels
cn[class %in% junklab, group_fine := "non-plankton/other"]

# ---- chain-forming flag from YOLO eligibility in the ROI data ---------------
d <- fread(cfg$roi_file, select = c("class_name","cell_count"), showProgress = FALSE)
elig <- d[, .(frac_chain = mean(cell_count != -1, na.rm = TRUE), n = .N), by = class_name]
cn <- merge(cn, elig, by.x = "class", by.y = "class_name", all.x = TRUE)
cn[, chain_forming := !is.na(frac_chain) & frac_chain >= cfg$chain_frac & group_fine == "diatom"]

# ---- coarse group used by the composition analysis (05) ----------------------
cn[, group := fcase(
  group_fine == "diatom" & chain_forming, "chain diatom",
  group_fine == "diatom",                 "diatom (other)",
  group_fine == "dinoflagellate",         "dinoflagellate",
  group_fine == "ciliate",                "ciliate",
  group_fine == "cryptophyte",            "cryptophyte",
  group_fine == "cyanobacteria",          "cyanobacteria",
  group_fine %in% c("chlorophyte","haptophyte","chrysophyte","euglenophyte","dictyochophyte","other stramenopile"),
                                          "other flagellate/nano",
  group_fine %in% c("heterotrophic flagellate","metazoan"), "heterotroph/other",
  default = "other")]

# ---- trophy from delivery taxonomy ------------------------------------------
cn[, trophy := fcase(trophic_type == "AU","autotroph", trophic_type == "MX","mixotroph",
                     trophic_type == "HT","heterotroph", trophic_type == "NS","unspecified",
                     default = "unknown")]

# ---- clean display name for figures -----------------------------------------
# species_flag_taxon_name, else class_clean ("Genus a Genus b" -> "Genus a/b"), else class
collapse_genus <- function(x) sub("^(\\S+) (\\S+) \\1 (\\S+)$", "\\1 \\2/\\3", x)
cn[, display_name := fifelse(!is.na(species_flag_taxon_name) & species_flag_taxon_name != "",
                             species_flag_taxon_name,
                     fifelse(!is.na(class_clean) & class_clean != "", collapse_genus(class_clean),
                             gsub("_", " ", class)))]

# ---- species flag on genus-level classes -------------------------------------
# Append spp./sp. where the display name lacks it: classes with an sflag but no composed
# name (e.g. Uroglenopsis_sp) and single-genus classes with no sflag (e.g. Eucampia).
# Multi-genus, cf. and species-level names are left as they are.
cn[, sflag_disp := fcase(toupper(sflag) == "SPP", "spp.",
                         toupper(sflag) == "SP",  "sp.",
                         rank == "Genus" & (is.na(sflag) | sflag == "") &
                           grepl("^[A-Z][a-z]+$", display_name), "spp.",
                         default = "")]
cn[sflag_disp != "" & !endsWith(display_name, sflag_disp),
   display_name := paste(display_name, sflag_disp)]

out <- cn[, .(class, display_name, scientificname, aphia_id, worms_class = wclass, base_group,
              group_fine, group, chain_forming, trophy, frac_chain_roi = round(frac_chain, 3))][order(group, class)]
fwrite(out, cfg$out_csv)

cat(sprintf("\nfunctional-group map -> %s (%d classes)\n", cfg$out_csv, nrow(out)))
cat("\n-- coarse group counts --\n"); print(out[, .N, by = group][order(-N)])
cat("\n-- fine group counts --\n");   print(out[, .N, by = group_fine][order(-N)])
cat(sprintf("\nchain-forming diatom classes: %d\n", out[chain_forming == TRUE, .N]))
print(out[chain_forming == TRUE, .(class, frac_chain_roi)][order(-frac_chain_roi)][1:12])
if (out[group_fine == "unassigned", .N] > 0) {
  cat("\nUnassigned classes, which need a manual taxonomy check:\n"); print(out[group_fine == "unassigned", .(class, scientificname, worms_class)])
}
