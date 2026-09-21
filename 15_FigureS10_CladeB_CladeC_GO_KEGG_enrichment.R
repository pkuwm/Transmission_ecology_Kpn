# Figure S10: GO and KEGG enrichment for Clade B and Clade C
rm(list = ls())
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(openxlsx)

# Settings
input_file <- file.choose()
output_dir <- file.path(dirname(input_file), "Enrichment_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clades <- c("CladeB", "CladeC")

clean_gene <- function(x) {
    x <- trimws(as.character(x))
    unique(x[!is.na(x) & x != ""])
}

gene_sets <- setNames(
    lapply(clades, function(clade) {
        clean_gene(read_excel(input_file, sheet = clade)$Gene)
    }),
    clades
)

# KEGG background
kegg <- read_excel(input_file, sheet = "KEGG背景库") %>%
    transmute(
        ID = trimws(as.character(Pathway)),
        Gene = trimws(as.character(Gene)),
        Description = as.character(Description)
    ) %>%
    filter(!is.na(ID), ID != "", !is.na(Gene), Gene != "") %>%
    distinct()

# GO background: the first two columns are Gene and GO_ID
go <- read_excel(input_file, sheet = "GO背景库")
go <- go[, 1:2]
names(go) <- c("Gene", "GO_ID")

go <- go %>%
    mutate(across(everything(), as.character)) %>%
    separate_rows(GO_ID, sep = "[,;]") %>%
    mutate(across(everything(), trimws)) %>%
    filter(
        !is.na(Gene), Gene != "",
        !is.na(GO_ID), GO_ID != ""
    ) %>%
    distinct()

valid_go <- intersect(
    unique(go$GO_ID),
    AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
)

if (length(valid_go) == 0) {
    stop("No valid GO identifiers found.")
}

go_annotation <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = valid_go,
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY", "DEFINITION")
) %>%
    transmute(
        GO_ID = GOID,
        Description = TERM,
        ONTOLOGY,
        Definition = DEFINITION
    )

go <- go %>%
    inner_join(go_annotation, by = "GO_ID") %>%
    rename(ID = GO_ID)

write.xlsx(
    go,
    file.path(output_dir, "GO_background.xlsx"),
    overwrite = TRUE
)

# Enrichment using genes represented in each annotation background
run_enrichment <- function(genes, background) {
    fit <- clusterProfiler::enricher(
        gene = genes,
        universe = unique(background$Gene),
        TERM2GENE = background %>%
            dplyr::select(ID, Gene) %>%
            distinct(),
        TERM2NAME = background %>%
            dplyr::select(ID, Description) %>%
            distinct(),
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500
    )
    
    result <- as.data.frame(fit)
    
    if (nrow(result) == 0) {
        return(result)
    }
    
    result %>%
        mutate(
            RichFactor = Count / as.numeric(sub("/.*", "", BgRatio))
        )
}

# Combined GO/KEGG bubble plot
plot_enrichment <- function(data, clade) {
    if (nrow(data) == 0) {
        message(clade, ": no terms passed the enrichment filters.")
        return(invisible(NULL))
    }
    
    ontology_colors <- c(
        BP = "#47706C",
        CC = "#ECB1CB",
        MF = "#C7C37F",
        KEGG = "#907287"
    )
    
    data <- data %>%
        mutate(
            ONTOLOGY = factor(ONTOLOGY, levels = names(ontology_colors)),
            Score = -log10(pmax(p.adjust, .Machine$double.xmin)),
            Term_key = paste(ONTOLOGY, ID, sep = ":")
        ) %>%
        arrange(ONTOLOGY, p.adjust) %>%
        mutate(Term_key = factor(Term_key, levels = rev(Term_key)))
    
    labels <- setNames(
        paste0(
            "<span style='color:",
            ontology_colors[as.character(data$ONTOLOGY)],
            "'>", data$Description, "</span>"
        ),
        as.character(data$Term_key)
    )
    
    p <- ggplot(data, aes(x = RichFactor, y = Term_key)) +
        geom_point(
            aes(size = Count, fill = Score),
            shape = 21,
            colour = "black",
            stroke = 0.35
        ) +
        scale_fill_gradientn(
            colours = c("#BAD2E5", "#95B0C4", "#6F8BA3"),
            name = "-log10(adjusted P)"
        ) +
        scale_size_continuous(
            range = c(3, 12),
            name = "Gene count"
        ) +
        scale_y_discrete(labels = labels) +
        scale_x_continuous(
            limits = c(0, max(data$RichFactor) * 1.2),
            expand = expansion(mult = 0)
        ) +
        labs(
            title = paste(
                sub("Clade", "Clade ", clade),
                "enrichment"
            ),
            x = "Rich factor",
            y = NULL
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.y = ggtext::element_markdown(size = 9),
            panel.grid = element_line(
                colour = "grey90",
                linetype = "dashed"
            )
        )
    
    ggsave(
        file.path(output_dir, paste0(clade, "_GO_KEGG_enrichment.pdf")),
        p,
        width = 8,
        height = max(5, 2 + 0.25 * nrow(data)),
        limitsize = FALSE,
        bg = "white",
        useDingbats = FALSE
    )
    
    invisible(p)
}

# Analyze both clades
results <- list()

for (clade in clades) {
    go_result <- run_enrichment(gene_sets[[clade]], go)
    kegg_result <- run_enrichment(gene_sets[[clade]], kegg)
    
    if (nrow(go_result) > 0) {
        go_result <- go_result %>%
            left_join(
                go %>% dplyr::select(ID, ONTOLOGY) %>% distinct(),
                by = "ID"
            )
    }
    
    if (nrow(kegg_result) > 0) {
        kegg_result$ONTOLOGY <- "KEGG"
    }
    
    results[[paste0(clade, "_GO")]] <- go_result
    results[[paste0(clade, "_KEGG")]] <- kegg_result
    
    plot_enrichment(bind_rows(go_result, kegg_result), clade)
}

write.xlsx(
    results,
    file.path(output_dir, "GO_KEGG_enrichment_results.xlsx"),
    overwrite = TRUE
)

message("Results saved in: ", output_dir)
