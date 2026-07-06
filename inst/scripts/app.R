library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(readxl)
library(specmine)
library(mclust)
library(ggplot2)

options(shiny.maxRequestSize = 100 * 1024^2)  # Allow uploads up to 100 MB

header_link <- function(inputId, label, icon_name = NULL) {
  icon_part <- if (!is.null(icon_name)) icon(icon_name) else NULL
  tags$li(
    class = "dropdown",
    actionLink(
      inputId = inputId,
      label = tagList(icon_part, label),
      style = "color:#f5f5f5; font-weight:600; padding:15px 12px; display:block; text-decoration:none;"
    )
  )
}

extract_pca_loadings <- function(pca_res, dims = 2) {
  if (!is.null(pca_res$rotation)) return(as.data.frame(pca_res$rotation[, seq_len(dims), drop = FALSE]))
  if (!is.null(pca_res$loadings)) return(as.data.frame(pca_res$loadings[, seq_len(dims), drop = FALSE]))
  if (!is.null(pca_res$loadings_matrix)) return(as.data.frame(pca_res$loadings_matrix[, seq_len(dims), drop = FALSE]))
  stop("Could not extract PCA loadings.")
}

build_loadings_plot <- function(loadings_df, pc = 1, top_n = 20) {
  pcname <- paste0("PC", pc)
  v <- loadings_df[[pcname]]
  ord <- order(abs(v), decreasing = TRUE)
  top_idx <- head(ord, top_n)
  
  df <- data.frame(
    Feature = rownames(loadings_df)[top_idx],
    Loading = v[top_idx],
    stringsAsFactors = FALSE
  )
  
  plotly::plot_ly(
    df,
    x = ~reorder(Feature, abs(Loading)),
    y = ~Loading,
    type = "bar"
  ) %>%
    plotly::layout(
      xaxis = list(title = "Feature"),
      yaxis = list(title = pcname)
    )
}

format_feature_labels <- function(features) {
  features <- as.character(features)
  if (all(grepl("^\\d+(\\.\\d+)?$", features))) {
    paste0(features, " cm-1")
  } else {
    features
  }
}

extract_pca_loadings <- function(pca_res, dims = 2) {
  if (!is.null(pca_res$rotation)) return(as.data.frame(pca_res$rotation[, seq_len(dims), drop = FALSE]))
  if (!is.null(pca_res$loadings)) return(as.data.frame(pca_res$loadings[, seq_len(dims), drop = FALSE]))
  if (!is.null(pca_res$loadings_matrix)) return(as.data.frame(pca_res$loadings_matrix[, seq_len(dims), drop = FALSE]))
  stop("Could not extract PCA loadings.")
}

top_features_table <- function(loadings_df, pcs = 1:3, top_n = 10) {
  pcs <- pcs[pcs <= ncol(loadings_df)]
  res <- lapply(pcs, function(pc) {
    v <- loadings_df[[pc]]
    ord <- order(abs(v), decreasing = TRUE)
    idx <- head(ord, top_n)
    data.frame(
      PC = paste0("PC", pc),
      Feature = rownames(loadings_df)[idx],
      FeatureLabel = format_feature_labels(rownames(loadings_df)[idx]),
      Loading = v[idx],
      AbsLoading = abs(v[idx]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

build_loadings_plot <- function(loadings_df, pc = 1, top_n = 20) {
  pcname <- paste0("PC", pc)
  v <- loadings_df[[pcname]]
  ord <- order(abs(v), decreasing = TRUE)
  idx <- head(ord, top_n)
  feats <- rownames(loadings_df)[idx]
  df <- data.frame(
    Feature = format_feature_labels(feats),
    RawFeature = feats,
    Loading = v[idx],
    stringsAsFactors = FALSE
  )
  plotly::plot_ly(
    df,
    x = ~reorder(Feature, abs(Loading)),
    y = ~Loading,
    type = "bar"
  ) %>%
    plotly::layout(
      xaxis = list(title = "Feature"),
      yaxis = list(title = pcname)
    )
}

build_loading_heatmap <- function(loadings_df, pcs = 1:3) {
  pcs <- pcs[pcs <= ncol(loadings_df)]
  mat <- as.matrix(loadings_df[, pcs, drop = FALSE])
  rownames(mat) <- format_feature_labels(rownames(mat))
  plotly::plot_ly(
    x = colnames(mat),
    y = rownames(mat),
    z = mat,
    type = "heatmap",
    colorscale = "RdBu",
    zmid = 0
  )
}

pca_interpretation_text <- function(loadings_df, pc = 1, top_n = 3) {
  pcname <- paste0("PC", pc)
  v <- loadings_df[[pcname]]
  ord <- order(abs(v), decreasing = TRUE)
  idx <- head(ord, top_n)
  feats <- format_feature_labels(rownames(loadings_df)[idx])
  vals <- round(v[idx], 3)
  paste0(
    pcname, " is mainly driven by ",
    paste(paste0(feats, " (", vals, ")"), collapse = ", "),
    "."
  )
}

build_loading_heatmap <- function(loadings_df, pcs = 1:3) {
  pcs <- pcs[pcs <= ncol(loadings_df)]
  mat <- as.matrix(loadings_df[, pcs, drop = FALSE])
  plotly::plot_ly(
    x = colnames(mat),
    y = rownames(mat),
    z = mat,
    type = "heatmap",
    colorscale = "RdBu",
    zmid = 0
  )
}

top_features_table <- function(loadings_df, pcs = 1:3, top_n = 10) {
  pcs <- pcs[pcs <= ncol(loadings_df)]
  res <- lapply(pcs, function(pc) {
    v <- loadings_df[[pc]]
    ord <- order(abs(v), decreasing = TRUE)
    idx <- head(ord, top_n)
    data.frame(
      PC = pc,
      Feature = rownames(loadings_df)[idx],
      Loading = v[idx],
      AbsLoading = abs(v[idx]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

build_hca_heatmap <- function(dat, scale_rows = TRUE) {
  x <- as.matrix(dat)
  if (scale_rows) x <- t(scale(t(x)))
  hc_rows <- hclust(dist(x))
  hc_cols <- hclust(dist(t(x)))
  list(matrix = x, row_hc = hc_rows, col_hc = hc_cols)
}

full_button <- function(id, label, icon_name = NULL, class = "btn-web") {
  actionButton(
    inputId = id,
    label = label,
    icon = if (!is.null(icon_name)) icon(icon_name) else NULL,
    class = class,
    width = "100%"
  )
}

safe_ggplotly <- function(p) {
  if (inherits(p, "plotly")) return(p)
  if (inherits(p, "gg") || inherits(p, "ggplot")) return(plotly::ggplotly(p))
  if (is.list(p)) {
    if (!is.null(p$plotly) && inherits(p$plotly, "plotly")) return(p$plotly)
    if (!is.null(p$plot) && inherits(p$plot, "plotly")) return(p$plot)
    if (!is.null(p$plot) && (inherits(p$plot, "gg") || inherits(p$plot, "ggplot"))) return(plotly::ggplotly(p$plot))
    if (!is.null(p$gg) && (inherits(p$gg, "gg") || inherits(p$gg, "ggplot"))) return(plotly::ggplotly(p$gg))
    if (length(p) >= 1 && inherits(p[[1]], "plotly")) return(p[[1]])
    if (length(p) >= 1 && (inherits(p[[1]], "gg") || inherits(p[[1]], "ggplot"))) return(plotly::ggplotly(p[[1]]))
  }
  stop("No applicable plot method for object returned by plotting function.")
}

read_uploaded_table <- function(fileinfo) {
  if (is.null(fileinfo)) return(NULL)
  
  ext <- tolower(tools::file_ext(fileinfo$name))
  
  if (ext == "csv") {
    df <- read.csv(fileinfo$datapath, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(fileinfo$datapath, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(readxl::read_excel(fileinfo$datapath))
  } else {
    stop("Unsupported file format. Use csv, tsv, txt, xlsx or xls.")
  }
  
  if (ncol(df) < 2) stop("The uploaded file must contain at least two columns.")
  
  first_col <- df[[1]]
  if (all(!is.na(first_col)) && anyDuplicated(first_col) == 0) {
    rownames(df) <- as.character(first_col)
    df <- df[, -1, drop = FALSE]
  }
  
  if (ncol(df) == 0) stop("No data columns remain after processing the first column as row names.")
  
  for (j in seq_len(ncol(df))) {
    if (is.character(df[[j]])) {
      suppressWarnings(num_col <- as.numeric(df[[j]]))
      if (!all(is.na(num_col))) df[[j]] <- num_col
    }
  }
  
  as.data.frame(df, check.names = FALSE)
}

make_default_metadata <- function(dat) {
  data.frame(
    Sample = colnames(dat),
    row.names = colnames(dat),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

create_dataset_safe <- function(dat, meta, selected_type = NULL) {
  dat <- as.data.frame(dat, check.names = FALSE, stringsAsFactors = FALSE)
  
  for (j in seq_len(ncol(dat))) {
    if (!is.numeric(dat[[j]])) {
      suppressWarnings(num_col <- as.numeric(dat[[j]]))
      if (!all(is.na(num_col))) {
        dat[[j]] <- num_col
      }
    }
  }
  
  ds <- tryCatch(
    create_dataset(
      data = as.matrix(dat),
      metadata = meta,
      type = selected_type
    ),
    error = function(e) NULL
  )
  
  if (is.null(ds)) {
    ds <- create_dataset(
      data = as.matrix(dat),
      metadata = meta
    )
  }
  
  ds
}

extract_pca_scores <- function(pca_res, dims = 2) {
  if (!is.null(pca_res$x)) return(as.data.frame(pca_res$x[, seq_len(dims), drop = FALSE]))
  if (!is.null(pca_res$scores)) return(as.data.frame(pca_res$scores[, seq_len(dims), drop = FALSE]))
  stop("Could not extract PCA scores from PCA result.")
}

get_embedding_matrix <- function(method, result_obj, dims = 2) {
  if (method == "PCA") {
    return(extract_pca_scores(result_obj, dims = dims))
  }
  if (!is.null(result_obj$embedding)) {
    return(as.data.frame(result_obj$embedding[, seq_len(dims), drop = FALSE]))
  }
  if (!is.null(result_obj$Y)) {
    return(as.data.frame(result_obj$Y[, seq_len(dims), drop = FALSE]))
  }
  stop(paste("Could not extract embedding coordinates for", method))
}

build_pca_plotly_3d <- function(pca_res, meta, class_col = NULL, show_labels = FALSE) {
  scores <- extract_pca_scores(pca_res, dims = 3)
  colnames(scores) <- c("PC1", "PC2", "PC3")
  
  meta2 <- meta[rownames(scores), , drop = FALSE]
  
  color_values <- if (!is.null(class_col) && class_col %in% colnames(meta2)) {
    as.factor(meta2[[class_col]])
  } else {
    as.factor(rownames(scores))
  }
  
  label_values <- rownames(scores)
  
  plot_ly(
    data = scores,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    type = "scatter3d",
    mode = if (isTRUE(show_labels)) "markers+text" else "markers",
    color = color_values,
    text = label_values,
    hoverinfo = "text",
    marker = list(size = 5)
  ) %>%
    layout(
      scene = list(
        xaxis = list(title = "PC1"),
        yaxis = list(title = "PC2"),
        zaxis = list(title = "PC3")
      ),
      legend = list(title = list(text = ifelse(is.null(class_col), "Group", class_col)))
    )
}

build_cluster_plot_pca_2d <- function(pca_res, cluster_vec, show_labels = FALSE) {
  scores <- extract_pca_scores(pca_res, dims = 2)
  colnames(scores) <- c("PC1", "PC2")
  scores$Cluster <- as.factor(cluster_vec)
  scores$Sample <- rownames(scores)
  
  plot_ly(
    data = scores,
    x = ~PC1,
    y = ~PC2,
    type = "scatter",
    mode = if (isTRUE(show_labels)) "markers+text" else "markers",
    color = ~Cluster,
    text = ~Sample,
    hoverinfo = "text",
    marker = list(size = 8)
  ) %>%
    layout(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2")
    )
}

replace_inf_with_na_df <- function(df) {
  x <- as.data.frame(df, check.names = FALSE)
  x[] <- lapply(x, function(col) {
    if (is.numeric(col)) {
      col[is.infinite(col)] <- NA
    }
    col
  })
  x
}

count_bad_values_df <- function(df) {
  x <- as.data.frame(df, check.names = FALSE)
  mat <- as.matrix(x)
  suppressWarnings(storage.mode(mat) <- "numeric")
  list(
    na = sum(is.na(mat)),
    nan = sum(is.nan(mat)),
    inf = sum(is.infinite(mat))
  )
}

safe_normalize <- function(ds, method = "sum") {
  tryCatch(
    normalize(ds, method = method),
    error = function(e) ds
  )
}

map_embedding_method <- function(x) {
  switch(
    x,
    "PCA" = "pca",
    "UMAP" = "umap",
    "t-SNE" = "tsne",
    "ICA" = "ica",
    tolower(x)
  )
}

get_builtin_dataset <- function(name) {
  if (name == "iris") {
    df <- datasets::iris
    meta <- data.frame(
      Species = df$Species,
      row.names = paste0("Sample_", seq_len(nrow(df))),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    dat <- t(as.matrix(df[, 1:4]))
    colnames(dat) <- rownames(meta)
    rownames(dat) <- colnames(df)[1:4]
    return(list(data = dat, metadata = meta, type = "metabolomics"))
  }
  
  if (name == "wine") {
    if (requireNamespace("rattle.data", quietly = TRUE)) {
      data("wine", package = "rattle.data")
      df <- get("wine", envir = asNamespace("rattle.data"))
    } else if (requireNamespace("rattle", quietly = TRUE)) {
      data("wine", package = "rattle")
      df <- get("wine", envir = asNamespace("rattle"))
    } else {
      stop("Package 'rattle.data' or 'rattle' is required for the built-in wine dataset.")
    }
    
    meta <- data.frame(
      Type = df$Type,
      row.names = paste0("Sample_", seq_len(nrow(df))),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    dat <- t(as.matrix(df[, setdiff(colnames(df), "Type"), drop = FALSE]))
    colnames(dat) <- rownames(meta)
    rownames(dat) <- setdiff(colnames(df), "Type")
    return(list(data = dat, metadata = meta, type = "metabolomics"))
  }
  
  if (name == "cachexia") {
    if (!requireNamespace("specmine.datasets", quietly = TRUE)) {
      stop("Package 'specmine.datasets' is required for the built-in cachexia dataset.")
    }
    
    data("cachexia", package = "specmine.datasets", envir = environment())
    obj <- get("cachexia", envir = environment())
    
    if (is.list(obj) && !is.null(obj$data) && !is.null(obj$metadata)) {
      dat <- obj$data
      meta <- obj$metadata
    } else {
      stop("The cachexia object does not contain both 'data' and 'metadata'.")
    }
    
    dat <- as.matrix(dat)
    meta <- as.data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)
    
    return(list(data = dat, metadata = meta, type = "metabolomics"))
  }
  
  if (name == "propolis") {
    if (!requireNamespace("specmine.datasets", quietly = TRUE)) {
      stop("Package 'specmine.datasets' is required for the built-in propolis dataset.")
    }
    
    data("propolis", package = "specmine.datasets", envir = environment())
    obj <- get("propolis", envir = environment())
    
    if (is.list(obj) && !is.null(obj$data) && !is.null(obj$metadata)) {
      dat <- obj$data
      meta <- obj$metadata
    } else {
      stop("The propolis object does not contain both 'data' and 'metadata'.")
    }
    
    dat <- as.matrix(dat)
    meta <- as.data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)
    
    return(list(data = dat, metadata = meta, type = "nmr-spectra"))
  }
  
  stop("Unknown built-in dataset.")
}

pcascoresplot2D <- function(dataset, pcares, column.class = NULL, labels = FALSE, ellipses = FALSE) {
  scores <- extract_pca_scores(pcares, dims = 2)
  colnames(scores) <- c("PC1", "PC2")
  scores$Sample <- rownames(scores)
  
  meta <- dataset$metadata
  if (!is.null(column.class) && !is.null(meta) && column.class %in% colnames(meta)) {
    scores$Group <- as.factor(meta[rownames(scores), column.class])
  } else {
    scores$Group <- factor("All")
  }
  
  p <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 3)
  
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = -0.8, size = 3, show.legend = FALSE)
  }
  
  if (isTRUE(ellipses)) {
    valid <- table(scores$Group)
    if (sum(valid >= 3) >= 2) {
      p <- p + ggplot2::stat_ellipse(type = "norm", linewidth = 0.8)
    }
  }
  
  p + ggplot2::theme_minimal() +
    ggplot2::labs(x = "PC1", y = "PC2", color = if (!is.null(column.class)) column.class else "Group")
}

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    titleWidth = 240,
    title = tags$div(
      style = "display:flex; align-items:center; gap:10px; font-weight:700;",
      tags$span(style = "font-size:22px;", "WebSpecmine")
    ),
    header_link("go_home", "Home", "home"),
    header_link("go_load_workspace", "Load Workspace", "folder-open"),
    header_link("go_preprocessing", "Preprocessing", "broom"),
    header_link("go_embeddings", "Embeddings", "chart-line"),
    header_link("go_clustering", "Clustering", "object-group"),
    header_link("go_compare", "Comparison", "columns"),
    header_link("go_results", "Results", "table")
  ),
  
  dashboardSidebar(
    width = 245,
    sidebarMenu(
      id = "tabs",
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Load Workspace", tabName = "load_workspace", icon = icon("folder-open")),
      menuItem("Preprocessing", tabName = "preprocessing", icon = icon("broom")),
      menuItem("Embeddings", tabName = "embeddings", icon = icon("project-diagram")),
      menuItem("Clustering", tabName = "clustering", icon = icon("object-group")),
      menuItem("Comparison", tabName = "comparison", icon = icon("columns")),
      menuItem("PCA / HCA", tabName = "pca_hca", icon = icon("project-diagram")),
      menuItem("Results", tabName = "results", icon = icon("table")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        body, .content-wrapper, .right-side {
          background: #eef2f7;
          font-family: Arial, Helvetica, sans-serif;
        }
        .skin-blue .main-header .logo {
          background: linear-gradient(90deg, #0a5b85, #0f8ca0);
          color: white;
          font-weight: 700;
          font-size: 24px;
          border-bottom: 0;
        }
        .skin-blue .main-header .logo:hover {
          background: linear-gradient(90deg, #0a5b85, #0f8ca0);
        }
        .skin-blue .main-header .navbar {
          background: linear-gradient(90deg, #0f8ca0, #10a7a7);
          border-bottom: 0;
        }
        .skin-blue .main-header .navbar .nav > li > a {
          color: #f5f5f5 !important;
          font-weight: 600;
        }
        .skin-blue .main-header .navbar .nav > li > a:hover,
        .skin-blue .main-header .navbar .nav > li > a:focus {
          background: rgba(255,255,255,0.10);
          color: white !important;
          text-decoration: none;
        }
        .skin-blue .main-sidebar {
          background-color: #050505;
        }
        .skin-blue .sidebar-menu > li > a {
          color: #e8e8e8;
          border-left: 3px solid transparent;
          font-size: 14px;
        }
        .skin-blue .sidebar-menu > li > a:hover {
          background: #111111;
          color: #ffffff;
          border-left-color: #17b7c9;
        }
        .skin-blue .sidebar-menu > li.active > a {
          background: #111111;
          color: #ffffff;
          border-left-color: #17b7c9;
          font-weight: 700;
        }
        .content {
          padding: 22px;
        }
        .webspecmine-banner {
          background: linear-gradient(90deg, #1198a6, #11aaa8);
          color: white;
          text-align: center;
          padding: 32px 20px 28px 20px;
          margin: -22px -22px 22px -22px;
          box-shadow: 0 2px 12px rgba(0,0,0,0.12);
        }
        .webspecmine-banner h1 {
          margin: 0;
          font-size: 34px;
          font-weight: 700;
        }
        .webspecmine-banner p {
          margin-top: 8px;
          margin-bottom: 0;
          font-size: 13px;
          text-transform: uppercase;
          letter-spacing: 1px;
        }
        .box {
          border-top: 0 !important;
          border-radius: 10px;
          box-shadow: 0 4px 14px rgba(0,0,0,0.08);
          overflow: hidden;
        }
        .box .box-header {
          padding-top: 14px;
          padding-bottom: 14px;
          background: #ffffff;
          border-bottom: 1px solid #edf1f5;
        }
        .box .box-title {
          font-weight: 700;
          color: #234;
        }
        .box .box-body {
          background: #ffffff;
        }
        .btn-web {
          background: linear-gradient(90deg, #2f80ed, #3c8df0);
          color: white !important;
          border: none;
          border-radius: 10px;
          padding: 10px 16px;
          font-weight: 700;
          box-shadow: 0 3px 10px rgba(60,141,240,0.18);
        }
        .btn-web:hover {
          background: linear-gradient(90deg, #2b76d8, #367fdd);
          color: white !important;
        }
        .btn-web-secondary {
          background: #12a6a6;
          color: white !important;
          border: none;
          border-radius: 10px;
          padding: 10px 16px;
          font-weight: 700;
          box-shadow: 0 3px 10px rgba(18,166,166,0.18);
        }
        .btn-web-secondary:hover {
          background: #0f9191;
          color: white !important;
        }
        .btn-web-danger {
          background: #d9534f;
          color: white !important;
          border: none;
          border-radius: 10px;
          padding: 10px 16px;
          font-weight: 700;
        }
        .btn-web-danger:hover {
          background: #c64541;
          color: white !important;
        }
        .small-note {
          color: #667085;
          font-size: 12px;
          line-height: 1.5;
        }
        .section-note {
          background: #f5fbfc;
          border-left: 4px solid #11aaa8;
          padding: 12px 14px;
          border-radius: 8px;
          color: #355;
          margin-top: 10px;
          margin-bottom: 4px;
          font-size: 13px;
        }
        .control-block {
          background: #f8fafc;
          border: 1px solid #edf2f7;
          border-radius: 10px;
          padding: 12px;
          margin-bottom: 14px;
        }
        .control-title {
          font-size: 13px;
          font-weight: 700;
          color: #344054;
          margin-bottom: 10px;
          text-transform: uppercase;
          letter-spacing: 0.3px;
        }
        .pre-status {
          background: #eef8ff;
          border: 1px solid #d7e9f8;
          border-radius: 10px;
          padding: 12px 14px;
          color: #245;
          margin-bottom: 10px;
        }
        .table .form-group {
          margin-bottom: 0;
        }
      "))
    ),
    
    tabItems(
      tabItem(
        tabName = "home",
        div(
          class = "webspecmine-banner",
          h1(icon("flask"), " WebSpecmine"),
          p("Metabolomics and Spectral Data Analysis and Mining")
        ),
        fluidRow(
          box(
            title = "Overview",
            width = 12,
            status = "primary",
            tags$p("This version keeps the main analysis workflow but now includes a stronger preprocessing step, an automatic cleaning shortcut, and a more polished interface.")
          )
        )
      ),
      
      tabItem(
        tabName = "load_workspace",
        fluidRow(
          box(
            title = "Built-in test dataset",
            width = 4,
            status = "primary",
            selectInput("builtin_dataset", "Choose built-in dataset", choices = c("iris", "wine", "cachexia", "propolis")),
            full_button("load_builtin", "Load built-in dataset", "database"),
            tags$div(class = "section-note", "Use this option for fast testing without uploading external files.")
          ),
          box(
            title = "Load Dataset",
            width = 4,
            status = "primary",
            fileInput("dataset_file", "Choose dataset", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            selectInput(
              "data_type",
              "Dataset type",
              choices = c("metabolomics", "spectra", "raman-spectra", "nmr-spectra", "ms-spectra"),
              selected = "metabolomics"
            ),
            full_button("load_data", "Load dataset", "upload"),
            tags$div(class = "section-note", "Rows should represent variables and columns should represent samples.")
          ),
          box(
            title = "Load Metadata",
            width = 4,
            status = "primary",
            fileInput("metadata_file", "Choose metadata (optional)", accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            tags$div(class = "section-note", "Metadata row names should match sample names.")
          )
        ),
        fluidRow(
          box(
            title = "Dataset Preview",
            width = 12,
            status = "primary",
            DTOutput("preview_table")
          )
        )
      ),
      
      tabItem(
        tabName = "preprocessing",
        fluidRow(
          valueBoxOutput("na_box", width = 4),
          valueBoxOutput("nan_box", width = 4),
          valueBoxOutput("inf_box", width = 4)
        ),
        fluidRow(
          box(
            title = "Auto clean",
            width = 4,
            status = "primary",
            div(
              class = "control-block",
              div(class = "control-title", "Automatic workflow"),
              checkboxInput("auto_replace_inf", "Convert Inf / -Inf to NA", TRUE),
              checkboxInput("auto_impute_median", "Impute missing values with median", TRUE),
              checkboxInput("auto_center", "Apply mean centering", FALSE),
              checkboxInput("auto_scale", "Apply scaling", TRUE),
              checkboxInput("auto_normalize", "Apply sum normalization", FALSE),
              full_button("auto_clean_btn", "Auto clean dataset", "magic", class = "btn-web-secondary")
            ),
            tags$div(class = "section-note", "Recommended shortcut before embeddings: clean invalid values, impute missing values, and optionally scale/center the dataset.")
          ),
          box(
            title = "Manual cleaning",
            width = 4,
            status = "primary",
            div(
              class = "control-block",
              div(class = "control-title", "Missing values and invalid values"),
              full_button("replace_inf_btn", "Convert Inf / -Inf to NA", "exchange-alt"),
              br(), br(),
              full_button("remove_var_na_btn", "Remove variables with NAs", "filter"),
              br(), br(),
              full_button("remove_sample_na_btn", "Remove samples with NAs", "user-times", class = "btn-web-danger")
            )
          ),
          box(
            title = "Imputation and transforms",
            width = 4,
            status = "primary",
            div(
              class = "control-block",
              div(class = "control-title", "Imputation"),
              selectInput(
                "impute_method",
                "Imputation method",
                choices = c("Mean" = "mean", "Median" = "median", "Fixed value" = "value", "kNN" = "knn")
              ),
              numericInput("impute_value", "Fixed value", value = 0),
              numericInput("impute_k", "k for kNN", value = 5, min = 1),
              full_button("impute_btn", "Apply imputation", "syringe")
            ),
            div(
              class = "control-block",
              div(class = "control-title", "Transformations"),
              checkboxInput("do_log", "Log transform", FALSE),
              checkboxInput("do_center", "Mean centering", FALSE),
              checkboxInput("do_scale", "Scaling", FALSE),
              checkboxInput("do_normalize", "Normalize by sum", FALSE),
              full_button("transform_btn", "Apply transformations", "sliders-h", class = "btn-web-secondary")
            )
          )
        ),
        fluidRow(
          box(
            title = "Preprocessing status",
            width = 12,
            status = "primary",
            div(class = "pre-status", verbatimTextOutput("preprocess_status"))
          )
        ),
        fluidRow(
          box(
            title = "Preprocessed dataset preview",
            width = 12,
            status = "primary",
            DTOutput("preprocess_preview_table")
          )
        )
      ),
      
      tabItem(
        tabName = "embeddings",
        fluidRow(
          box(
            title = "Embedding Controls",
            width = 3,
            status = "primary",
            selectInput("embedding_method", "Method", choices = c("PCA", "UMAP", "t-SNE", "ICA")),
            selectInput("embedding_color", "Colour by", choices = "Sample"),
            radioButtons(
              "embedding_dims",
              "Dimensions",
              choices = c("2D" = 2, "3D" = 3),
              selected = 2,
              inline = TRUE
            ),
            checkboxInput("embedding_scale", "Scale before analysis", FALSE),
            checkboxInput("embedding_center", "Center before analysis", TRUE),
            numericInput("umap_neighbors", "UMAP n_neighbors", value = 15, min = 2),
            numericInput("umap_min_dist", "UMAP min_dist", value = 0.1, min = 0, max = 1, step = 0.05),
            textInput("umap_metric", "UMAP metric", value = "euclidean"),
            numericInput("tsne_perplexity", "t-SNE perplexity", value = 30, min = 2),
            numericInput("ica_maxit", "ICA max iterations", value = 200, min = 50),
            checkboxInput("embedding_labels", "Show labels", FALSE),
            checkboxInput("pca_ellipses", "PCA ellipses", FALSE),
            full_button("run_embedding", "Run embedding", "play")
          ),
          box(
            title = "Embedding Plot",
            width = 9,
            status = "primary",
            plotlyOutput("embedding_plot", height = 520)
          )
        ),
        fluidRow(
          box(
            title = "Embedding Status",
            width = 12,
            status = "primary",
            verbatimTextOutput("embedding_status")
          )
        )
      ),
      
      tabItem(
        tabName = "clustering",
        fluidRow(
          box(
            title = "Clustering Controls",
            width = 3,
            status = "primary",
            selectInput("clustering_method", "Method", choices = c("DBSCAN", "HDBSCAN", "GMM", "kmeans", "hc")),
            selectInput("cluster_plot_method", "Plot coordinates", choices = c("PCA", "UMAP", "t-SNE", "ICA")),
            selectInput("cluster_dims", "Plot dimensions", choices = c("2D", "3D"), selected = "2D"),
            numericInput("dbscan_eps", "DBSCAN eps", value = 0.5, min = 0.01, step = 0.05),
            numericInput("dbscan_minPts", "DBSCAN/HDBSCAN minPts", value = 5, min = 2),
            numericInput("gmm_clusters", "GMM number of clusters", value = 3, min = 1),
            numericInput("kmeans_clusters", "k-means number of clusters", value = 3, min = 1),
            selectInput("hc_distance", "HC distance", choices = c("euclidean", "manhattan", "pearson", "spearman")),
            selectInput("hc_linkage", "HC linkage", choices = c("complete", "average", "single", "ward.D", "ward.D2")),
            checkboxInput("clustering_scale", "Scale before clustering", FALSE),
            full_button("run_clustering", "Run clustering", "project-diagram")
          ),
          box(
            title = "Clustering Plot",
            width = 9,
            status = "primary",
            plotlyOutput("clustering_plot", height = 520)
          )
        ),
        fluidRow(
          box(
            title = "Clustering Status",
            width = 12,
            status = "primary",
            verbatimTextOutput("clustering_status")
          )
        )
      ),
      
      tabItem(
        tabName = "pca_hca",
        fluidRow(
          box(
            title = "PCA / HCA controls",
            width = 3,
            status = "primary",
            numericInput("pca_top_n", "Top features per PC", value = 15, min = 5),
            numericInput("pca_heatmap_pcs", "Number of PCs in heatmap", value = 3, min = 2),
            checkboxInput("hca_scale_rows", "Scale rows for heatmap", TRUE),
            numericInput("hca_clusters", "HCA cut clusters", value = 3, min = 2),
            full_button("run_pca_hca", "Run PCA/HCA", "project-diagram")
          ),
          box(
            title = "PCA Loadings",
            width = 9,
            status = "primary",
            plotlyOutput("pca_loadings_plot", height = 350)
          )
        ),
        fluidRow(
          box(
            title = "Interpretation",
            width = 12,
            status = "primary",
            verbatimTextOutput("pca_interpretation")
          )
        ),
        fluidRow(
          box(
            title = "Loading Heatmap",
            width = 6,
            status = "primary",
            plotlyOutput("pca_loading_heatmap", height = 450)
          ),
          box(
            title = "HCA Dendrogram",
            width = 6,
            status = "primary",
            plotOutput("hca_dendrogram", height = 450)
          )
        ),
        fluidRow(
          box(
            title = "HCA Heatmap",
            width = 6,
            status = "primary",
            plotlyOutput("hca_heatmap", height = 450)
          ),
          box(
            title = "Top Bands per PC",
            width = 6,
            status = "primary",
            DTOutput("top_features_table")
          )
        )
      ),
      
      tabItem(
        tabName = "comparison",
        fluidRow(
          box(
            title = "Comparison Controls",
            width = 3,
            status = "primary",
            checkboxGroupInput(
              "compare_embeddings_select",
              "Embeddings to compare",
              choices = c("PCA", "UMAP", "t-SNE", "ICA"),
              selected = c("PCA", "UMAP", "t-SNE")
            ),
            numericInput("compare_neighbors", "Quality metric neighbours", value = 5, min = 2),
            full_button("run_compare_embeddings", "Compare embeddings", "chart-bar"),
            br(), br(),
            checkboxGroupInput(
              "compare_clusterings_select",
              "Clusterings to compare",
              choices = c("DBSCAN", "HDBSCAN", "GMM", "kmeans", "hc"),
              selected = c("DBSCAN", "HDBSCAN", "GMM")
            ),
            full_button("run_compare_clusterings", "Compare clusterings", "balance-scale")
          ),
          box(
            title = "Embedding Comparison",
            width = 9,
            status = "primary",
            DTOutput("embedding_comparison_table")
          )
        ),
        fluidRow(
          box(
            title = "Clustering Comparison",
            width = 12,
            status = "primary",
            DTOutput("clustering_comparison_table")
          )
        )
      ),
      
      tabItem(
        tabName = "results",
        fluidRow(
          valueBoxOutput("trust_box", width = 4),
          valueBoxOutput("cont_box", width = 4),
          valueBoxOutput("sil_box", width = 4)
        ),
        fluidRow(
          box(
            title = "Current Result Summary",
            width = 12,
            status = "primary",
            verbatimTextOutput("result_summary")
          )
        )
      ),
      
      tabItem(
        tabName = "help",
        fluidRow(
          box(
            title = "Help",
            width = 12,
            status = "primary",
            tags$p("1. Upload a csv/tsv/txt/xlsx file or load a built-in test dataset."),
            tags$p("2. Use Preprocessing to inspect NA / NaN / Inf values."),
            tags$p("3. Use Auto clean for a quick workflow, or the manual controls for more detailed preparation."),
            tags$p("4. Run PCA, UMAP, t-SNE or ICA."),
            tags$p("5. Run DBSCAN, HDBSCAN, GMM, k-means or hierarchical clustering."),
            tags$p("6. Compare methods in the Comparison tab.")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(input$go_home, { updateTabItems(session, "tabs", "home") })
  observeEvent(input$go_load_workspace, { updateTabItems(session, "tabs", "load_workspace") })
  observeEvent(input$go_preprocessing, { updateTabItems(session, "tabs", "preprocessing") })
  observeEvent(input$go_embeddings, { updateTabItems(session, "tabs", "embeddings") })
  observeEvent(input$go_clustering, { updateTabItems(session, "tabs", "clustering") })
  observeEvent(input$go_compare, { updateTabItems(session, "tabs", "comparison") })
  observeEvent(input$go_results, { updateTabItems(session, "tabs", "results") })
  
  rv <- reactiveValues(
    raw_data = NULL,
    metadata = NULL,
    dataset = NULL,
    dataset_type = NULL,
    embedding_results = list(),
    clustering_results = list(),
    embedding_quality = NULL,
    clustering_quality = NULL,
    current_source = NULL,
    preprocess_log = "No preprocessing applied yet."
  )
  
  rv$pca_hca_results <- reactiveValues()
  
  observeEvent(input$run_pca_hca, {
    req(rv$dataset)
    tryCatch({
      pca_res <- pca_analysis_dataset(rv$dataset, scale = TRUE, center = TRUE)
      rv$pca_hca_results$pca <- pca_res
      rv$pca_hca_results$loadings <- extract_pca_loadings(pca_res, dims = max(3, input$pca_heatmap_pcs))
      rv$pca_hca_results$hca <- hclust(dist(t(as.matrix(rv$dataset$data))))
      showNotification("PCA/HCA analysis completed.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("PCA/HCA error:", e$message), type = "error", duration = NULL)
    })
  })
  
  output$pca_interpretation <- renderText({
    req(rv$pca_hca_results$loadings)
    pcs <- seq_len(min(3, input$pca_heatmap_pcs))
    paste(vapply(
      pcs,
      function(pc) pca_interpretation_text(rv$pca_hca_results$loadings, pc = pc, top_n = 3),
      character(1)
    ), collapse = "\n")
  })
  
  output$pca_loadings_plot <- renderPlotly({
    req(rv$pca_hca_results$loadings)
    build_loadings_plot(rv$pca_hca_results$loadings, pc = 1, top_n = input$pca_top_n)
  })
  
  output$pca_loading_heatmap <- renderPlotly({
    req(rv$pca_hca_results$loadings)
    build_loading_heatmap(rv$pca_hca_results$loadings, pcs = seq_len(input$pca_heatmap_pcs))
  })
  
  output$hca_dendrogram <- renderPlot({
    req(rv$pca_hca_results$hca)
    plot(as.dendrogram(rv$pca_hca_results$hca), main = "HCA Dendrogram", ylab = "Height")
  })
  
  output$hca_heatmap <- renderPlotly({
    req(rv$dataset)
    res <- build_hca_heatmap(rv$dataset$data, scale_rows = input$hca_scale_rows)
    plotly::plot_ly(
      x = colnames(res$matrix),
      y = rownames(res$matrix),
      z = res$matrix,
      type = "heatmap",
      colorscale = "RdBu",
      zmid = 0
    )
  })
  
  output$top_features_table <- renderDT({
    req(rv$pca_hca_results$loadings)
    datatable(
      top_features_table(rv$pca_hca_results$loadings, pcs = seq_len(input$pca_heatmap_pcs), top_n = input$pca_top_n),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  
  rebuild_dataset <- function() {
    req(rv$raw_data, rv$metadata)
    rv$dataset <- create_dataset_safe(
      dat = rv$raw_data,
      meta = rv$metadata,
      selected_type = rv$dataset_type
    )
    rv$embedding_results <- list()
    rv$clustering_results <- list()
    rv$embedding_quality <- NULL
    rv$clustering_quality <- NULL
  }
  
  reset_analysis_results <- function() {
    rv$embedding_results <- list()
    rv$clustering_results <- list()
    rv$embedding_quality <- NULL
    rv$clustering_quality <- NULL
  }
  
  output$preview_table <- renderDT({
    if (!is.null(rv$raw_data)) {
      datatable(
        head(as.data.frame(rv$raw_data, check.names = FALSE), 10),
        options = list(scrollX = TRUE, pageLength = 5),
        rownames = TRUE
      )
    } else {
      req(input$dataset_file)
      tmp <- tryCatch(
        read_uploaded_table(input$dataset_file),
        error = function(e) {
          showNotification(paste("Preview error:", e$message), type = "error", duration = NULL)
          NULL
        }
      )
      req(tmp)
      datatable(
        head(tmp, 10),
        options = list(scrollX = TRUE, pageLength = 5),
        rownames = TRUE
      )
    }
  })
  
  output$preprocess_preview_table <- renderDT({
    req(rv$raw_data)
    datatable(
      head(as.data.frame(rv$raw_data, check.names = FALSE), 10),
      options = list(scrollX = TRUE, pageLength = 5),
      rownames = TRUE
    )
  })
  
  observeEvent(input$load_builtin, {
    tryCatch({
      obj <- get_builtin_dataset(input$builtin_dataset)
      
      rv$raw_data <- as.data.frame(obj$data, check.names = FALSE)
      rv$metadata <- obj$metadata
      rv$dataset_type <- obj$type
      rv$dataset <- create_dataset_safe(
        dat = obj$data,
        meta = obj$metadata,
        selected_type = obj$type
      )
      reset_analysis_results()
      rv$current_source <- paste("Built-in:", input$builtin_dataset)
      rv$preprocess_log <- "Dataset loaded. No preprocessing applied yet."
      
      updateSelectInput(
        session,
        "embedding_color",
        choices = colnames(obj$metadata),
        selected = colnames(obj$metadata)[1]
      )
      
      showNotification(
        paste("Built-in dataset", input$builtin_dataset, "loaded successfully."),
        type = "message",
        duration = 5
      )
    }, error = function(e) {
      showNotification(
        paste("Error loading built-in dataset:", e$message),
        type = "error",
        duration = NULL
      )
    })
  })
  
  observeEvent(input$load_data, {
    req(input$dataset_file)
    
    tryCatch({
      dat <- read_uploaded_table(input$dataset_file)
      
      meta <- if (!is.null(input$metadata_file)) {
        read_uploaded_table(input$metadata_file)
      } else {
        make_default_metadata(dat)
      }
      
      if (!all(colnames(dat) %in% rownames(meta))) {
        meta <- make_default_metadata(dat)
        showNotification(
          "Metadata names did not match sample names. Default metadata was created.",
          type = "warning",
          duration = 8
        )
      } else {
        meta <- meta[colnames(dat), , drop = FALSE]
      }
      
      rv$raw_data <- dat
      rv$metadata <- meta
      rv$dataset_type <- input$data_type
      rv$dataset <- create_dataset_safe(
        dat = dat,
        meta = meta,
        selected_type = input$data_type
      )
      reset_analysis_results()
      rv$current_source <- paste("Uploaded:", input$dataset_file$name)
      rv$preprocess_log <- "Dataset loaded. No preprocessing applied yet."
      
      updateSelectInput(
        session,
        "embedding_color",
        choices = colnames(meta),
        selected = if (ncol(meta) > 0) colnames(meta)[1] else "Sample"
      )
      
      showNotification("Dataset loaded successfully.", type = "message", duration = 5)
      
    }, error = function(e) {
      rv$raw_data <- NULL
      rv$metadata <- NULL
      rv$dataset <- NULL
      rv$dataset_type <- NULL
      
      showNotification(
        paste("Error loading dataset:", e$message),
        type = "error",
        duration = NULL
      )
    })
  })
  
  observeEvent(input$replace_inf_btn, {
    req(rv$raw_data)
    tryCatch({
      rv$raw_data <- replace_inf_with_na_df(rv$raw_data)
      rebuild_dataset()
      rv$preprocess_log <- "Manual preprocessing: Inf and -Inf values converted to NA."
      showNotification("Inf and -Inf values converted to NA.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Preprocessing error:", e$message), type = "error", duration = NULL)
    })
  })
  
  observeEvent(input$remove_var_na_btn, {
    req(rv$dataset)
    tryCatch({
      rv$dataset <- remove_variables_by_nas(rv$dataset)
      rv$raw_data <- as.data.frame(rv$dataset$data, check.names = FALSE)
      rv$metadata <- as.data.frame(rv$dataset$metadata, check.names = FALSE)
      reset_analysis_results()
      rv$preprocess_log <- "Manual preprocessing: variables with missing values removed."
      showNotification("Variables with missing values removed.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Preprocessing error:", e$message), type = "error", duration = NULL)
    })
  })
  
  observeEvent(input$remove_sample_na_btn, {
    req(rv$dataset)
    tryCatch({
      rv$dataset <- remove_samples_by_nas(rv$dataset)
      rv$raw_data <- as.data.frame(rv$dataset$data, check.names = FALSE)
      rv$metadata <- as.data.frame(rv$dataset$metadata, check.names = FALSE)
      reset_analysis_results()
      
      updateSelectInput(
        session,
        "embedding_color",
        choices = colnames(rv$metadata),
        selected = if (ncol(rv$metadata) > 0) colnames(rv$metadata)[1] else character(0)
      )
      
      rv$preprocess_log <- "Manual preprocessing: samples with missing values removed."
      showNotification("Samples with missing values removed.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Preprocessing error:", e$message), type = "error", duration = NULL)
    })
  })
  
  observeEvent(input$impute_btn, {
    req(rv$dataset)
    tryCatch({
      ds <- rv$dataset
      
      if (input$impute_method == "mean") {
        ds <- impute_nas_mean(ds)
      } else if (input$impute_method == "median") {
        ds <- impute_nas_median(ds)
      } else if (input$impute_method == "value") {
        ds <- impute_nas_value(ds, value = input$impute_value)
      } else if (input$impute_method == "knn") {
        ds <- impute_nas_knn(ds, k = input$impute_k)
      }
      
      rv$dataset <- ds
      rv$raw_data <- as.data.frame(ds$data, check.names = FALSE)
      rv$metadata <- as.data.frame(ds$metadata, check.names = FALSE)
      reset_analysis_results()
      rv$preprocess_log <- paste("Manual preprocessing: imputation applied using", input$impute_method, ".")
      
      showNotification("Imputation applied successfully.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Imputation error:", e$message), type = "error", duration = NULL)
    })
  })
  
  observeEvent(input$transform_btn, {
    req(rv$dataset)
    tryCatch({
      ds <- rv$dataset
      
      if (isTRUE(input$do_log)) {
        ds <- log_transform(ds)
      }
      if (isTRUE(input$do_center)) {
        ds <- mean_centering(ds)
      }
      if (isTRUE(input$do_scale)) {
        ds <- scaling(ds)
      }
      if (isTRUE(input$do_normalize)) {
        ds <- safe_normalize(ds, method = "sum")
      }
      
      rv$dataset <- ds
      rv$raw_data <- as.data.frame(ds$data, check.names = FALSE)
      rv$metadata <- as.data.frame(ds$metadata, check.names = FALSE)
      reset_analysis_results()
      
      applied <- c(
        if (isTRUE(input$do_log)) "log transform" else NULL,
        if (isTRUE(input$do_center)) "mean centering" else NULL,
        if (isTRUE(input$do_scale)) "scaling" else NULL,
        if (isTRUE(input$do_normalize)) "sum normalization" else NULL
      )
      rv$preprocess_log <- paste("Manual preprocessing:", paste(applied, collapse = ", "), ".")
      
      showNotification("Transformations applied successfully.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Transformation error:", e$message), type = "error", duration = NULL)
    })
  })
  
  observeEvent(input$auto_clean_btn, {
    req(rv$dataset)
    tryCatch({
      ds <- rv$dataset
      steps <- c()
      
      if (isTRUE(input$auto_replace_inf)) {
        rv$raw_data <- replace_inf_with_na_df(rv$raw_data)
        ds <- create_dataset_safe(rv$raw_data, rv$metadata, rv$dataset_type)
        steps <- c(steps, "Inf/-Inf converted to NA")
      }
      
      if (isTRUE(input$auto_impute_median)) {
        ds <- impute_nas_median(ds)
        steps <- c(steps, "missing values imputed with median")
      }
      
      if (isTRUE(input$auto_center)) {
        ds <- mean_centering(ds)
        steps <- c(steps, "mean centering")
      }
      
      if (isTRUE(input$auto_scale)) {
        ds <- scaling(ds)
        steps <- c(steps, "scaling")
      }
      
      if (isTRUE(input$auto_normalize)) {
        ds <- safe_normalize(ds, method = "sum")
        steps <- c(steps, "sum normalization")
      }
      
      rv$dataset <- ds
      rv$raw_data <- as.data.frame(ds$data, check.names = FALSE)
      rv$metadata <- as.data.frame(ds$metadata, check.names = FALSE)
      reset_analysis_results()
      
      if (length(steps) == 0) {
        rv$preprocess_log <- "Auto clean ran, but no step was selected."
      } else {
        rv$preprocess_log <- paste("Auto clean applied:", paste(steps, collapse = "; "), ".")
      }
      
      showNotification("Auto clean completed successfully.", type = "message", duration = 5)
    }, error = function(e) {
      showNotification(paste("Auto clean error:", e$message), type = "error", duration = NULL)
    })
  })
  
  output$preprocess_status <- renderText({
    req(rv$raw_data)
    cnt <- count_bad_values_df(rv$raw_data)
    paste(
      "Source:", if (is.null(rv$current_source)) "None" else rv$current_source,
      "\nPreprocessing log:", rv$preprocess_log,
      "\nCurrent NA count:", cnt$na,
      "\nCurrent NaN count:", cnt$nan,
      "\nCurrent Inf count:", cnt$inf,
      "\nSamples:", ncol(rv$raw_data),
      "\nVariables:", nrow(rv$raw_data)
    )
  })
  
  output$na_box <- renderValueBox({
    req(rv$raw_data)
    cnt <- count_bad_values_df(rv$raw_data)
    valueBox(as.character(cnt$na), "NA values", icon = icon("question-circle"), color = "yellow")
  })
  
  output$nan_box <- renderValueBox({
    req(rv$raw_data)
    cnt <- count_bad_values_df(rv$raw_data)
    valueBox(as.character(cnt$nan), "NaN values", icon = icon("exclamation-circle"), color = "orange")
  })
  
  output$inf_box <- renderValueBox({
    req(rv$raw_data)
    cnt <- count_bad_values_df(rv$raw_data)
    valueBox(as.character(cnt$inf), "Inf values", icon = icon("infinity"), color = "red")
  })
  
  observeEvent(input$run_embedding, {
    req(rv$dataset)
    
    tryCatch({
      method <- input$embedding_method
      dims <- as.numeric(input$embedding_dims)
      
      res <- switch(
        method,
        "PCA" = pca_analysis_dataset(
          rv$dataset,
          scale = input$embedding_scale,
          center = input$embedding_center
        ),
        "UMAP" = umap_analysis_dataset(
          rv$dataset,
          n_components = dims,
          n_neighbors = input$umap_neighbors,
          min_dist = input$umap_min_dist,
          metric = input$umap_metric,
          scale = input$embedding_scale
        ),
        "t-SNE" = tsne_analysis_dataset(
          rv$dataset,
          n_components = dims,
          perplexity = input$tsne_perplexity,
          scale = input$embedding_scale
        ),
        "ICA" = ica_analysis_dataset(
          rv$dataset,
          n_components = dims,
          maxit = input$ica_maxit,
          scale = input$embedding_scale
        )
      )
      
      rv$embedding_results[[method]] <- res
      
      emb_for_metrics <- get_embedding_matrix(method, res, dims = min(2, dims))
      
      rv$embedding_quality <- embedding_quality_metrics(
        dataset = rv$dataset,
        embedding = emb_for_metrics,
        scale = input$embedding_scale,
        n_neighbors = 5
      )
      
      showNotification(paste(method, "completed."), type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(
        paste("Error running embedding:", e$message),
        type = "error",
        duration = NULL
      )
    })
  })
  
  output$embedding_plot <- renderPlotly({
    req(rv$dataset)
    req(rv$embedding_results[[input$embedding_method]])
    
    method <- input$embedding_method
    dims <- as.numeric(input$embedding_dims)
    class_col <- input$embedding_color
    
    build_embedding_plot_2d <- function(emb_res, meta, method, class_col = NULL, show_labels = FALSE) {
      emb <- get_embedding_matrix(method, emb_res, dims = 2)
      colnames(emb) <- c("X1", "X2")
      emb$Sample <- rownames(emb)
      
      if (!is.null(class_col) && class_col %in% colnames(meta)) {
        emb$Group <- as.factor(meta[rownames(emb), class_col])
      } else {
        emb$Group <- as.factor("All samples")
      }
      
      plotly::plot_ly(
        data = emb,
        x = ~X1,
        y = ~X2,
        type = "scatter",
        mode = if (isTRUE(show_labels)) "markers+text" else "markers",
        color = ~Group,
        text = ~Sample,
        hoverinfo = "text",
        marker = list(size = 8)
      ) %>%
        plotly::layout(
          xaxis = list(title = paste(method, "1")),
          yaxis = list(title = paste(method, "2"))
        )
    }
    
    p <- switch(
      method,
      "PCA" = {
        if (dims == 2) {
          build_embedding_plot_2d(
            rv$embedding_results[[method]],
            rv$metadata,
            "PCA",
            class_col,
            input$embedding_labels
          )
        } else {
          build_pca_plotly_3d(
            pca_res = rv$embedding_results[[method]],
            meta = rv$metadata,
            class_col = class_col,
            show_labels = input$embedding_labels
          )
        }
      },
      "UMAP" = {
        if (dims == 2) {
          build_embedding_plot_2d(
            rv$embedding_results[[method]],
            rv$metadata,
            "UMAP",
            class_col,
            input$embedding_labels
          )
        } else {
          umap_scoresplot3D(rv$dataset, rv$embedding_results[[method]], column.class = class_col)
        }
      },
      "t-SNE" = {
        if (dims == 2) {
          build_embedding_plot_2d(
            rv$embedding_results[[method]],
            rv$metadata,
            "t-SNE",
            class_col,
            input$embedding_labels
          )
        } else {
          tsne_scoresplot3D(rv$dataset, rv$embedding_results[[method]], column.class = class_col)
        }
      },
      "ICA" = {
        if (dims == 2) {
          build_embedding_plot_2d(
            rv$embedding_results[[method]],
            rv$metadata,
            "ICA",
            class_col,
            input$embedding_labels
          )
        } else {
          ica_scoresplot3D(rv$dataset, rv$embedding_results[[method]], column.class = class_col)
        }
      }
    )
    
    safe_ggplotly(p)
  })
  
  output$embedding_status <- renderText({
    req(rv$dataset)
    paste(
      "Source:", if (is.null(rv$current_source)) "None" else rv$current_source,
      "\nLoaded dataset:", !is.null(rv$dataset),
      "\nAvailable embeddings:", if (length(names(rv$embedding_results)) > 0) paste(names(rv$embedding_results), collapse = ", ") else "None",
      "\nSelected method:", input$embedding_method,
      "\nDimensions:", as.numeric(input$embedding_dims)
    )
  })
  
  observeEvent(input$run_clustering, {
    req(rv$dataset)
    
    tryCatch({
      method <- input$clustering_method
      
      res <- switch(
        method,
        "DBSCAN" = dbscan_analysis_dataset(
          rv$dataset,
          eps = input$dbscan_eps,
          minPts = input$dbscan_minPts,
          scale = input$clustering_scale
        ),
        "HDBSCAN" = hdbscan_analysis_dataset(
          rv$dataset,
          minPts = input$dbscan_minPts,
          scale = input$clustering_scale
        ),
        "GMM" = gmm_analysis_dataset(
          rv$dataset,
          num.clusters = input$gmm_clusters,
          scale = input$clustering_scale
        ),
        "kmeans" = clustering(
          dataset = rv$dataset,
          method = "kmeans",
          num.clusters = input$kmeans_clusters
        ),
        "hc" = clustering(
          dataset = rv$dataset,
          method = "hc",
          distance = input$hc_distance,
          clustMethod = input$hc_linkage
        )
      )
      
      rv$clustering_results[[method]] <- res
      
      if (method %in% c("DBSCAN", "HDBSCAN", "GMM") && !is.null(res$cluster)) {
        uniq_cl <- unique(res$cluster)
        uniq_cl <- uniq_cl[!is.na(uniq_cl)]
        uniq_cl_no_noise <- setdiff(uniq_cl, -1)
        
        if (length(uniq_cl_no_noise) >= 2) {
          rv$clustering_quality <- cluster_quality_metrics(
            dataset = rv$dataset,
            clusters = res$cluster,
            scale = input$clustering_scale,
            remove.noise = TRUE
          )
        } else {
          rv$clustering_quality <- NULL
          showNotification(
            "At least 2 clusters are required to compute clustering quality metrics.",
            type = "warning",
            duration = 6
          )
        }
      } else {
        rv$clustering_quality <- NULL
      }
      
      showNotification(paste(method, "completed."), type = "message", duration = 5)
      
    }, error = function(e) {
      showNotification(
        paste("Error running clustering:", e$message),
        type = "error",
        duration = NULL
      )
    })
  })
  
  output$clustering_plot <- renderPlotly({
    req(rv$dataset)
    req(rv$clustering_results[[input$clustering_method]])
    
    cl_method <- input$clustering_method
    coord_method <- input$cluster_plot_method
    use3d <- identical(input$cluster_dims, "3D")
    coord_method_specmine <- map_embedding_method(coord_method)
    
    get_cluster_vec <- function(res, n = NULL) {
      if (is.null(res)) stop("Clustering result is NULL.")
      
      candidates <- c("cluster", "clusters", "clustering", "memberships", "membership")
      for (nm in candidates) {
        if (!is.null(res[[nm]])) return(res[[nm]])
      }
      
      if (inherits(res, "hclust")) {
        k <- n
        if (is.null(k)) k <- 2
        return(cutree(res, k = k))
      }
      
      if (is.atomic(res) && is.numeric(res)) return(res)
      if (is.factor(res)) return(as.integer(res))
      
      stop(paste(
        "No cluster vector found in clustering result. Available fields:",
        paste(names(res), collapse = ", ")
      ))
    }
    
    plot_cluster_from_pca <- function(pca_res, cluster_vec, show_labels = FALSE) {
      scores <- extract_pca_scores(pca_res, dims = 2)
      scores <- as.data.frame(scores)
      
      if (nrow(scores) == 0) stop("PCA scores are empty.")
      if (ncol(scores) < 2) stop("PCA scores do not contain 2 dimensions.")
      if (length(cluster_vec) != nrow(scores)) {
        stop(paste("Cluster vector length", length(cluster_vec), "does not match PCA rows", nrow(scores)))
      }
      
      colnames(scores)[1:2] <- c("PC1", "PC2")
      scores$Cluster <- as.factor(cluster_vec)
      scores$Sample <- rownames(scores)
      
      plotly::plot_ly(
        data = scores,
        x = ~PC1,
        y = ~PC2,
        type = "scatter",
        mode = if (isTRUE(show_labels)) "markers+text" else "markers",
        color = ~Cluster,
        text = ~Sample,
        hoverinfo = "text",
        marker = list(size = 8)
      ) %>%
        plotly::layout(
          xaxis = list(title = "PC1"),
          yaxis = list(title = "PC2")
        )
    }
    
    p <- switch(
      cl_method,
      "DBSCAN" = {
        emb_res <- rv$embedding_results[[coord_method]]
        req(emb_res)
        
        if (coord_method == "PCA") {
          validate(need(!use3d, "DBSCAN PCA plot is available only in 2D here."))
          plot_cluster_from_pca(
            pca_res = emb_res,
            cluster_vec = get_cluster_vec(rv$clustering_results[[cl_method]]),
            show_labels = input$embedding_labels
          )
        } else if (!use3d) {
          dbscan_plot2D(
            rv$dataset,
            rv$clustering_results[[cl_method]],
            method = coord_method_specmine,
            result.obj = emb_res,
            scale = input$clustering_scale
          )
        } else {
          dbscan_plot3D(
            rv$dataset,
            rv$clustering_results[[cl_method]],
            method = coord_method_specmine,
            result.obj = emb_res
          )
        }
      },
      "HDBSCAN" = {
        emb_res <- rv$embedding_results[[coord_method]]
        req(emb_res)
        
        if (coord_method == "PCA") {
          validate(need(!use3d, "HDBSCAN PCA plot is available only in 2D here."))
          plot_cluster_from_pca(
            pca_res = emb_res,
            cluster_vec = get_cluster_vec(rv$clustering_results[[cl_method]]),
            show_labels = input$embedding_labels
          )
        } else if (!use3d) {
          hdbscan_plot2D(
            rv$dataset,
            rv$clustering_results[[cl_method]],
            method = coord_method_specmine,
            result.obj = emb_res,
            scale = input$clustering_scale
          )
        } else {
          hdbscan_plot3D(
            rv$dataset,
            rv$clustering_results[[cl_method]],
            method = coord_method_specmine,
            result.obj = emb_res
          )
        }
      },
      "GMM" = {
        emb_res <- rv$embedding_results[[coord_method]]
        req(emb_res)
        
        if (coord_method == "PCA") {
          validate(need(!use3d, "GMM PCA plot is available only in 2D here."))
          plot_cluster_from_pca(
            pca_res = emb_res,
            cluster_vec = get_cluster_vec(rv$clustering_results[[cl_method]]),
            show_labels = input$embedding_labels
          )
        } else if (!use3d) {
          gmm_plot2D(
            rv$dataset,
            rv$clustering_results[[cl_method]],
            method = coord_method_specmine,
            result.obj = emb_res,
            scale = input$clustering_scale
          )
        } else {
          gmm_plot3D(
            rv$dataset,
            rv$clustering_results[[cl_method]],
            method = coord_method_specmine,
            result.obj = emb_res
          )
        }
      },
      "kmeans" = {
        req(rv$embedding_results[["PCA"]])
        validate(need(!use3d, "k-means plot is available only in 2D here."))
        pca_kmeans_plot2D(
          rv$dataset,
          rv$embedding_results[["PCA"]],
          num.clusters = input$kmeans_clusters,
          kmeans.result = rv$clustering_results[[cl_method]]
        )
      },
      "hc" = {
        req(rv$embedding_results[["PCA"]])
        validate(need(!use3d, "Hierarchical clustering plot is available only in 2D here."))
        cl_vec <- get_cluster_vec(rv$clustering_results[[cl_method]], n = input$kmeans_clusters)
        plot_cluster_from_pca(
          pca_res = rv$embedding_results[["PCA"]],
          cluster_vec = cl_vec,
          show_labels = input$embedding_labels
        )
      }
    )
    
    safe_ggplotly(p)
  })
  
  output$clustering_status <- renderText({
    req(rv$dataset)
    paste(
      "Available clusterings:", if (length(names(rv$clustering_results)) > 0) paste(names(rv$clustering_results), collapse = ", ") else "None",
      "\nSelected clustering:", input$clustering_method,
      "\nCoordinate method:", input$cluster_plot_method,
      "\nPlot dimensions:", input$cluster_dims
    )
  })
  
  embedding_comparison_res <- eventReactive(input$run_compare_embeddings, {
    req(rv$dataset)
    
    selected <- input$compare_embeddings_select
    validate(need(length(selected) > 0, "Select at least one embedding."))
    
    emb_list <- list()
    
    if ("PCA" %in% selected) {
      if (is.null(rv$embedding_results[["PCA"]])) {
        rv$embedding_results[["PCA"]] <- pca_analysis_dataset(rv$dataset)
      }
      emb_list[["PCA"]] <- get_embedding_matrix("PCA", rv$embedding_results[["PCA"]], dims = 2)
    }
    
    if ("UMAP" %in% selected) {
      if (is.null(rv$embedding_results[["UMAP"]])) {
        rv$embedding_results[["UMAP"]] <- umap_analysis_dataset(rv$dataset, n_components = 2)
      }
      emb_list[["UMAP"]] <- get_embedding_matrix("UMAP", rv$embedding_results[["UMAP"]], dims = 2)
    }
    
    if ("t-SNE" %in% selected) {
      if (is.null(rv$embedding_results[["t-SNE"]])) {
        rv$embedding_results[["t-SNE"]] <- tsne_analysis_dataset(rv$dataset, n_components = 2)
      }
      emb_list[["t-SNE"]] <- get_embedding_matrix("t-SNE", rv$embedding_results[["t-SNE"]], dims = 2)
    }
    
    if ("ICA" %in% selected) {
      if (is.null(rv$embedding_results[["ICA"]])) {
        rv$embedding_results[["ICA"]] <- ica_analysis_dataset(rv$dataset, n_components = 2)
      }
      emb_list[["ICA"]] <- get_embedding_matrix("ICA", rv$embedding_results[["ICA"]], dims = 2)
    }
    
    compare_embeddings(
      dataset = rv$dataset,
      embeddings.list = emb_list,
      scale = FALSE,
      n_neighbors = input$compare_neighbors
    )
  })
  
  output$embedding_comparison_table <- renderDT({
    req(embedding_comparison_res())
    datatable(
      embedding_comparison_res(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  clustering_comparison_res <- eventReactive(input$run_compare_clusterings, {
    req(rv$dataset)
    
    selected <- input$compare_clusterings_select
    validate(need(length(selected) > 0, "Select at least one clustering."))
    
    cl_list <- list()
    
    if ("DBSCAN" %in% selected) {
      if (is.null(rv$clustering_results[["DBSCAN"]])) {
        rv$clustering_results[["DBSCAN"]] <- dbscan_analysis_dataset(rv$dataset)
      }
      cl_list[["DBSCAN"]] <- rv$clustering_results[["DBSCAN"]]$cluster
    }
    
    if ("HDBSCAN" %in% selected) {
      if (is.null(rv$clustering_results[["HDBSCAN"]])) {
        rv$clustering_results[["HDBSCAN"]] <- hdbscan_analysis_dataset(rv$dataset)
      }
      cl_list[["HDBSCAN"]] <- rv$clustering_results[["HDBSCAN"]]$cluster
    }
    
    if ("GMM" %in% selected) {
      if (is.null(rv$clustering_results[["GMM"]])) {
        rv$clustering_results[["GMM"]] <- gmm_analysis_dataset(rv$dataset, num.clusters = 3)
      }
      cl_list[["GMM"]] <- rv$clustering_results[["GMM"]]$cluster
    }
    
    if ("kmeans" %in% selected) {
      if (is.null(rv$clustering_results[["kmeans"]])) {
        rv$clustering_results[["kmeans"]] <- clustering(rv$dataset, method = "kmeans", num.clusters = 3)
      }
      cl_list[["kmeans"]] <- rv$clustering_results[["kmeans"]]$cluster
    }
    
    compare_clusterings(
      dataset = rv$dataset,
      clusterings.list = cl_list,
      scale = FALSE,
      remove.noise = TRUE
    )
  })
  
  output$clustering_comparison_table <- renderDT({
    req(clustering_comparison_res())
    datatable(
      clustering_comparison_res(),
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$trust_box <- renderValueBox({
    val <- if (!is.null(rv$embedding_quality) && !is.null(rv$embedding_quality$trustworthiness)) {
      round(rv$embedding_quality$trustworthiness, 3)
    } else {
      NA
    }
    valueBox(as.character(val), "Trustworthiness", icon = icon("check-circle"), color = "aqua")
  })
  
  output$cont_box <- renderValueBox({
    val <- if (!is.null(rv$embedding_quality) && !is.null(rv$embedding_quality$continuity)) {
      round(rv$embedding_quality$continuity, 3)
    } else {
      NA
    }
    valueBox(as.character(val), "Continuity", icon = icon("sync"), color = "blue")
  })
  
  output$sil_box <- renderValueBox({
    val <- if (!is.null(rv$clustering_quality) && !is.null(rv$clustering_quality$silhouette)) {
      round(rv$clustering_quality$silhouette, 3)
    } else {
      NA
    }
    valueBox(as.character(val), "Silhouette", icon = icon("chart-bar"), color = "teal")
  })
  
  output$result_summary <- renderText({
    emb_names <- names(rv$embedding_results)
    cl_names <- names(rv$clustering_results)
    cnt <- if (!is.null(rv$raw_data)) count_bad_values_df(rv$raw_data) else list(na = NA, nan = NA, inf = NA)
    
    paste(
      "Source:", if (is.null(rv$current_source)) "None" else rv$current_source,
      "\nEmbeddings run:", if (length(emb_names) == 0) "None" else paste(emb_names, collapse = ", "),
      "\nClusterings run:", if (length(cl_names) == 0) "None" else paste(cl_names, collapse = ", "),
      "\nEmbedding quality available:", !is.null(rv$embedding_quality),
      "\nClustering quality available:", !is.null(rv$clustering_quality),
      "\nCurrent NA count:", cnt$na,
      "\nCurrent NaN count:", cnt$nan,
      "\nCurrent Inf count:", cnt$inf,
      "\nPreprocessing log:", rv$preprocess_log
    )
  })
}

shinyApp(ui, server)