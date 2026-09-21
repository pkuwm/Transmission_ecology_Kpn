# SNP cutpoint estimation, clusters and putative transmission events
rm(list = ls())
library(readxl)
library(dplyr)
library(writexl)
library(cutpointr)
library(igraph)

# Settings
snp_dir <- choose.dir(caption = "Select the SNP matrix folder")
if (is.na(snp_dir)) stop("No folder selected.")

message("Select the isolate metadata Excel file.")
info_file <- file.choose()

threshold <- 21
output_dir <- file.path(
    snp_dir,
    paste0("SNP_results_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read matrices and retain each unordered pair once
snp_files <- list.files(
    snp_dir,
    pattern = "SWHOsnp\\.csv$",
    full.names = TRUE
)

if (length(snp_files) == 0) stop("No SNP matrices found.")

snp_pairs <- bind_rows(lapply(snp_files, function(file) {
    m <- as.matrix(read.csv(
        file,
        row.names = 1,
        check.names = FALSE
    ))
    
    rownames(m) <- trimws(rownames(m))
    colnames(m) <- trimws(colnames(m))
    
    if (
        nrow(m) != ncol(m) ||
        !setequal(rownames(m), colnames(m))
    ) {
        stop("Matrix row and column names do not match: ", basename(file))
    }
    
    m <- m[, rownames(m), drop = FALSE]
    storage.mode(m) <- "numeric"
    
    if (
        anyNA(m) ||
        any(!is.finite(m)) ||
        any(m < 0) ||
        any(diag(m) != 0) ||
        !isTRUE(all.equal(m, t(m)))
    ) {
        stop("Invalid or asymmetric SNP matrix: ", basename(file))
    }
    
    index <- which(upper.tri(m), arr.ind = TRUE)
    
    data.frame(
        Isolate1 = rownames(m)[index[, 1]],
        Isolate2 = colnames(m)[index[, 2]],
        SNP_distance = m[index],
        ST_type = rep(
            sub("SWHOsnp\\.csv$", "", basename(file)),
            nrow(index)
        )
    )
}))

# Read metadata
info <- read_excel(info_file, sheet = 1) %>%
    select(StrainID, Province, Location, Niche) %>%
    mutate(across(everything(), ~ trimws(as.character(.x))))

if (anyDuplicated(info$StrainID)) {
    stop("Duplicate StrainID values in metadata.")
}

info1 <- info %>%
    rename(
        Isolate1 = StrainID,
        Province1 = Province,
        Location1 = Location,
        Niche1 = Niche
    )

info2 <- info %>%
    rename(
        Isolate2 = StrainID,
        Province2 = Province,
        Location2 = Location,
        Niche2 = Niche
    )

pairs <- snp_pairs %>%
    left_join(info1, by = "Isolate1") %>%
    left_join(info2, by = "Isolate2")

metadata_columns <- c(
    "Province1", "Location1", "Niche1",
    "Province2", "Location2", "Niche2"
)

if (any(vapply(
    pairs[metadata_columns],
    function(x) any(is.na(x) | x == ""),
    logical(1)
))) {
    stop("Missing metadata or unmatched isolate IDs.")
}

pairs <- pairs %>%
    mutate(
        NicheTransmission = if_else(
            Province1 == Province2 &
                Location1 == Location2 &
                Niche1 == Niche2,
            "same",
            "different"
        )
    )

# Estimate the pooled cutpoint
screening_pairs <- pairs %>%
    filter(SNP_distance <= 100)

cutpoint_result <- cutpointr(
    screening_pairs,
    x = SNP_distance,
    class = NicheTransmission,
    pos_class = "different",
    method = maximize_metric,
    metric = sum_sens_spec,
    boot_runs = 0
)

print(cutpoint_result)
plot(cutpoint_result)

pdf(
    file.path(output_dir, "01_cutpoint_plot.pdf"),
    width = 8,
    height = 5
)
tryCatch(plot(cutpoint_result), finally = dev.off())


# Identify clusters using the fixed SNP threshold
linked_pairs <- pairs %>%
    filter(SNP_distance <= threshold)

if (nrow(linked_pairs) == 0) {
    stop("No linked pairs found.")
}

g <- graph_from_data_frame(
    linked_pairs %>% select(Isolate1, Isolate2),
    directed = FALSE
)

membership <- components(g)$membership

membership_table <- data.frame(
    StrainID = names(membership),
    Cluster = as.integer(membership)
) %>%
    left_join(info, by = "StrainID") %>%
    arrange(Cluster, StrainID)

cluster_pairs <- linked_pairs %>%
    left_join(
        membership_table %>%
            select(Isolate1 = StrainID, Cluster),
        by = "Isolate1"
    ) %>%
    select(Cluster, everything()) %>%
    arrange(Cluster, SNP_distance)

write_xlsx(
    list(
        membership = membership_table,
        linked_pairs = cluster_pairs
    ),
    file.path(output_dir, "02_clusters_and_pairs.xlsx")
)

# Define one event per cluster and unordered institution pair
# Location identifiers must be unique across provinces
events <- cluster_pairs %>%
    mutate(
        Institution_A = pmin(Location1, Location2),
        Institution_B = pmax(Location1, Location2),
        Event_Type = if_else(
            Location1 == Location2,
            "Within-Institution",
            "Cross-Institution"
        )
    ) %>%
    distinct(Cluster, Event_Type, Institution_A, Institution_B) %>%
    arrange(Cluster, Event_Type, Institution_A, Institution_B) %>%
    mutate(
        Event_ID = sprintf("E%03d", row_number()),
        Transmission_Location = paste(
            Institution_A, Institution_B, sep = " - "
        )
    ) %>%
    select(
        Event_ID, Cluster, Event_Type,
        Institution_A, Institution_B, Transmission_Location
    )

write_xlsx(
    events,
    file.path(output_dir, "03_transmission_events.xlsx")
)

message(
    "Linked pairs: ", nrow(cluster_pairs),
    "; Clusters: ", n_distinct(membership_table$Cluster),
    "; Events: ", nrow(events)
)

message("Results saved in: ", output_dir)
