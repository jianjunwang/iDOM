
#' @title Optimized specialization index H2' of DOM-microbe associations
#' @description This function calculates the network-level specialization index H2' for full, negative, and positive DOM-microbe bipartite networks.
#' This optimized version keeps the original `H2()` function unchanged and adds optional sample-level parallel execution.
#' @param mol.data Sample/DOM matrix with samples in rows and DOM molecules in columns (compositional or abundance data).
#' @param micro.data Sample/microbe matrix with samples in rows and microbial taxa in columns (compositional or abundance data).
#' @param occu.rate The threshold for retaining bacterial taxa or DOM molecules observed in a specified proportion of the total samples. Default: 0.3.
#' @param sparcc.R An optional precomputed SparCC correlation matrix. If `NULL`, it will be computed internally. Default: `NULL`.
#' @param sparcc.R.threshold The absolute SparCC threshold used to retain DOM-microbe correlations in the bipartite network. Default: 0.3.
#' @param N Number of null models to generate. Default: 100.
#' @param Null.model Null model type. Can be given as an integer or name: 1/"r2dtable", 2/"swap.web", 3/"vaznull", 4/"shuffle.web". Default: `"shuffle.web"`.
#' @param network.type Character vector giving the network types to analyse. Available options are `"Full"`, `"Positive"`, and `"Negative"`. Default is all three.
#' @param ncores Number of CPU cores used for sample-level parallel computing. Default is 3, which uses serial execution.
#' @param verbose Logical; if `TRUE`, emits status messages. Default: `interactive()`. When enabled, `H2.fast()` reports preprocessing and SparCC status messages together with network-analysis progress updates.
#' @return A data frame containing `Index`, `Observed`, `Standardised`, `p.value`, `Network.type`, and `Sample`.
#' @details
#' The goal of `H2.fast()` is to accelerate the existing workflow without changing the core H2' calculation, null-model strategy, or output structure used by `H2()`.
#' If `sparcc.R` is not supplied, the suggested package `SpiecEasi` is required.
#' When `ncores > 1`, sample-level calculations are distributed with a cross-platform PSOCK cluster via `parallel::parLapply()`.
#' @references A Hu, M Choi, A J Tanentzap, J Liu, K-S Jang, J T Lennon, Y Liu, J Soininen, X Lu, Y Zhang, J Shen, and J Wang. 2022. 
#' Ecological networks of dissolved organic matter and microorganisms under global change. 
#' *Nature Communications*. [https://www.nature.com/articles/s41467-022-31251-1](https://www.nature.com/articles/s41467-022-31251-1)
#' 
#' F Meng, A Hu, K-S Jang, and J Wang. 2025. 
#' iDOM: Statistical analysis of dissolved organic matter characterized by high-resolution mass spectrometry. 
#' *mLife*. [https://onlinelibrary.wiley.com/doi/10.1002/mlf2.70002](https://onlinelibrary.wiley.com/doi/10.1002/mlf2.70002)
#' 
#' @examples
#' \donttest{
#' data(mol.data)
#' data(micro.data)
#'
#' mol_ex <- mol.data[, 1:60]
#' mic_ex <- micro.data[, 1:60]
#'
#' norm.mol <- vegan::decostand(mol_ex, method = "total")
#' norm.mic <- vegan::decostand(mic_ex, method = "total")
#' mm <- cbind(
#'   as.matrix(round(100000 * norm.mol, 0)),
#'   as.matrix(round(100000 * norm.mic, 0))
#' )
#' keep <- colSums(mm > 0) >= nrow(mm) * 0.3
#' mm_keep <- mm[, keep, drop = FALSE]
#' if (requireNamespace("SpiecEasi", quietly = TRUE)) {
#'   sparcc_ex <- SpiecEasi::sparcc(mm_keep)[["Cor"]]
#'   colnames(sparcc_ex) <- rownames(sparcc_ex) <- colnames(mm_keep)
#' } else {
#'   sparcc_ex <- stats::cor(mm_keep)
#' }
#'
#' res <- H2.fast(
#'   mol.data = mol_ex,
#'   micro.data = mic_ex,
#'   sparcc.R = sparcc_ex,
#'   N = 5,
#'   ncores = 1,
#'   network.type = "Full",
#'   verbose = FALSE
#' )
#' head(res)
#' }
#' @seealso
#' \code{\link{H2}}, \code{\link[vegan]{decostand}},
#' \code{\link[bipartite]{nullmodel}}, \code{\link[bipartite]{networklevel}}
#' @export

H2.fast <- function(mol.data,
                    micro.data,
                    occu.rate = 0.3,
                    sparcc.R = NULL,
                    sparcc.R.threshold = 0.3,
                    N = 100,
                    Null.model = "shuffle.web",
                    network.type = c("Full", "Positive", "Negative"),
                    ncores = 3,
                    verbose = interactive()) {
  
  if (length(N) != 1 || is.na(N) || N < 1 || N != as.integer(N)) {
    stop("N must be a positive integer.")
  }
  if (length(occu.rate) != 1 || is.na(occu.rate) || occu.rate < 0 || occu.rate > 1) {
    stop("occu.rate must be between 0 and 1.")
  }
  network.type <- unique(as.character(network.type))
  valid_network_types <- c("Full", "Positive", "Negative")
  if (length(network.type) < 1 || any(!network.type %in% valid_network_types)) {
    stop("network.type must contain one or more of: Full, Positive, Negative.")
  }
  verbose <- isTRUE(verbose)
  progress_state <- list(
    enabled = verbose,
    start_time = proc.time()[["elapsed"]],
    last_label = NULL
  )
  
  update_progress <- function(done, label) {
    if (isTRUE(progress_state$enabled)) {
      .h2_progress_update(
        progress_state,
        completed_units = done,
        total_units = total_units,
        label = label
      )
    }
  }
  
  emit_status <- function(...) {
    if (isTRUE(progress_state$enabled)) {
      message(...)
    }
  }
  
  make_result <- function(sample, class_name, observed = NA_real_, standardised = NA_real_, p_value = NA_real_) {
    data.frame(
      Index = "H2",
      Observed = observed,
      Standardised = standardised,
      p.value = p_value,
      Network.type = class_name,
      Sample = sample
    )
  }
  
  # 1. Standardize the two abundance tables and keep features that pass the occurrence filter across all samples.
  norm.mol.data <- vegan::decostand(mol.data, method = "total")
  norm.micro.data <- vegan::decostand(micro.data, method = "total")
  norm.mol.data.round <- as.matrix(round(100000 * norm.mol.data, 0))
  norm.micro.data.round <- as.matrix(round(100000 * norm.micro.data, 0))
  if (!identical(rownames(norm.mol.data.round), rownames(norm.micro.data.round))) {
    stop("Row names of mol.data and micro.data do not match.")
  }
  mol.micro.data <- cbind(norm.mol.data.round, norm.micro.data.round)
  keep_mask <- colSums(mol.micro.data > 0) >= nrow(mol.micro.data) * occu.rate
  mol.micro.data.keep <- mol.micro.data[, keep_mask, drop = FALSE]
  if (ncol(mol.micro.data.keep) == 0) {
    stop("No features were retained after applying occu.rate. Lower occu.rate and try again.")
  }
  emit_status("Preprocessing complete (", ncol(mol.micro.data.keep), " retained features)")

  # 2. Obtain the SparCC correlation matrix, either by using the supplied matrix or by computing it from the filtered abundance table.
  compute_sparcc <- is.null(sparcc.R)
  if (compute_sparcc) {
    emit_status("Running SparCC correlation estimation on ", ncol(mol.micro.data.keep), " retained features")
  }
  if (compute_sparcc) {
    check_suggested("SpiecEasi", "H2.fast")
    sparcc.out <- SpiecEasi::sparcc(data = mol.micro.data.keep)
    sparcc.R <- sparcc.out$Cor
    colnames(sparcc.R) <- rownames(sparcc.R) <- colnames(mol.micro.data.keep)
  } else {
    if (!identical(colnames(mol.micro.data.keep), rownames(sparcc.R))) {
      stop("Column names of mol.micro.data.keep do not match row names of sparcc.R.")
    }
    if (!identical(colnames(mol.micro.data.keep), colnames(sparcc.R))) {
      stop("Column names of mol.micro.data.keep do not match column names of sparcc.R.")
    }
    sparcc.R <- as.matrix(sparcc.R)
  }

  # 3. Build the DOM-microbe interaction matrix and precompute, for each sample, which DOM molecules and microbes are present after filtering.
  sparcc.R.thresh <- ifelse(abs(sparcc.R) >= sparcc.R.threshold, sparcc.R, 0)
  dom_keep_names <- intersect(colnames(mol.data), colnames(mol.micro.data.keep))
  micro_keep_names <- intersect(colnames(micro.data), colnames(mol.micro.data.keep))
  sparcc.R.thresh.DOM.Micro <- sparcc.R.thresh[dom_keep_names, micro_keep_names, drop = FALSE]
  DOM.Micro.int <- round(100000 * sparcc.R.thresh.DOM.Micro, 0)
  
  sample_names <- rownames(mol.micro.data.keep)
  sample_present <- lapply(seq_len(nrow(mol.micro.data.keep)), function(jj) {
    list(
      dom = dom_keep_names[mol.micro.data.keep[jj, dom_keep_names, drop = TRUE] > 0],
      micro = micro_keep_names[mol.micro.data.keep[jj, micro_keep_names, drop = TRUE] > 0]
    )
  })

  # 4. Prepare the selected network types and the optional parallel backend.
  class_mats <- list(
    Full = abs(DOM.Micro.int),
    Positive = ifelse(DOM.Micro.int > 0, DOM.Micro.int, 0),
    Negative = ifelse(DOM.Micro.int < 0, abs(DOM.Micro.int), 0)
  )
  class_mats <- class_mats[network.type]
  
  emit_status(if (compute_sparcc) "SparCC correlation estimation complete" else "Using precomputed SparCC matrix")
  
  # 5. Prepare the execution scaffold for the selected network types. This includes the chosen network names, optional parallel backend, and progress bookkeeping used during the sample-by-sample calculations below.
  class_names <- names(class_mats)
  parallel_info <- .h2_prepare_parallel(ncores, nrow(mol.micro.data.keep))
  total_samples <- parallel_info$total_samples
  sample_ids <- parallel_info$sample_ids
  sample_chunks <- parallel_info$sample_chunks
  cl <- parallel_info$cl
  total_units <- total_samples * length(class_mats)
  if (!is.null(cl)) {
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

  # 5. Analyse each requested network type separately.
  result_chunks <- vector("list", length(class_mats))
  completed_units <- 0L
  
  for (class_idx in seq_along(class_mats)) {
    class_name <- class_names[[class_idx]]
    adj.cor.int <- class_mats[[class_idx]]
    row_keep <- rowSums(adj.cor.int) > 0
    col_keep <- colSums(adj.cor.int) > 0
    row_keep_idx <- which(row_keep)
    col_keep_idx <- which(col_keep)

    if (!any(row_keep) || !any(col_keep)) {
      result_chunks[[class_idx]] <- do.call(
        rbind,
        lapply(sample_names, make_result, class_name = class_name)
      )
      completed_units <- completed_units + total_samples
      update_progress(completed_units, paste0(class_name, " network has no retained links; skipping"))
      next
    }

    # Compress the current network to the rows/columns that still have links.
    adj.cor.int.class <- adj.cor.int[row_keep, col_keep, drop = FALSE]

    # Generate null networks once for the current network type.
    nulls <- bipartite::nullmodel(adj.cor.int.class, N = N, method = Null.model)
    nulls <- lapply(nulls, function(x) {
      dimnames(x) <- dimnames(adj.cor.int.class)
      x
    })

    # For each sample, find which DOM molecules and microbes still remain in the current class-specific network, then convert them to row/column positions in the compact matrix.
    sample_dom_names_class <- lapply(sample_present, function(x) {
      intersect(x$dom, rownames(adj.cor.int.class))
    })
    sample_micro_names_class <- lapply(sample_present, function(x) {
      intersect(x$micro, colnames(adj.cor.int.class))
    })
    
    compute_one <- function(jj) {
      dom_names <- sample_dom_names_class[[jj]]
      micro_names <- sample_micro_names_class[[jj]]

      if (length(dom_names) == 0L || length(micro_names) == 0L) {
        return(make_result(sample_names[[jj]], class_name))
      }

      adj.cor.int.site <- adj.cor.int.class[dom_names, micro_names, drop = FALSE]
      index.obs <- as.numeric(bipartite::networklevel(
        adj.cor.int.site, weighted = TRUE, index = "H2"
      ))

      index.null <- vapply(nulls, function(null_mat) {
        as.numeric(bipartite::networklevel(
          null_mat[dom_names, micro_names, drop = FALSE], weighted = TRUE, index = "H2"
        ))
      }, numeric(1))

      index.rand.mean <- mean(index.null)
      index.rand.sd <- stats::sd(index.null)
      index.ses <- if (is.na(index.rand.sd) || index.rand.sd == 0) {
        NA_real_
      } else {
        (index.obs - index.rand.mean) / index.rand.sd
      }

      praw <- sum(index.null > index.obs) / length(index.null)
      index.obs.p <- ifelse(praw > 0.5, 1 - praw, praw)

      make_result(
        sample = sample_names[[jj]], class_name = class_name, observed = index.obs, standardised = index.ses, p_value = index.obs.p
      )
    }

    class_res <- vector("list", length(sample_ids))
    if (!is.null(cl)) {
      # In parallel mode, workers receive the current network, its null models, and the sample-to-index mappings for this network type.
      parallel::clusterExport(
        cl,
        varlist = c(
          "sample_names", "sample_dom_names_class", "sample_micro_names_class",
          "adj.cor.int.class", "nulls", "class_name",
          "make_result",
          "compute_one"
        ),
        envir = environment()
      )
      offset <- 0L
      for (chunk in sample_chunks) {
        chunk_res <- parallel::parLapply(
          cl, chunk, compute_one
        )
        class_res[(offset + 1L):(offset + length(chunk_res))] <- chunk_res
        offset <- offset + length(chunk_res)
        update_progress(
          completed_units + offset,
          paste0(class_name, " network: ", offset, "/", total_samples, " samples")
        )
      }
    } else {
      for (kk in seq_along(sample_ids)) {
        class_res[[kk]] <- compute_one(sample_ids[[kk]])
        update_progress(
          completed_units + kk,
          paste0(class_name, " network: ", kk, "/", total_samples, " samples")
        )
      }
    }
    result_chunks[[class_idx]] <- do.call(rbind, class_res)
    completed_units <- completed_units + total_samples
  }

  update_progress(total_units, "H2.fast complete")
  do.call(rbind, result_chunks)
}

#' @noRd
.h2_prepare_parallel <- function(ncores, total_samples) {
  ncores <- as.integer(ncores[[1]])
  if (is.na(ncores) || ncores < 1L) {
    ncores <- 1L
  }

  detected <- parallel::detectCores(logical = FALSE)
  if (is.na(detected) || detected < 1L) {
    detected <- ncores
  }

  ncores <- min(ncores, detected, total_samples)
  sample_ids <- seq_len(total_samples)
  chunk_size <- max(1L, ceiling(length(sample_ids) / max(1L, ncores * 2L)))
  sample_chunks <- split(sample_ids, ceiling(seq_along(sample_ids) / chunk_size))

  cl <- NULL
  if (ncores > 1L) {
    cl <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(cl, {
      loadNamespace("bipartite")
      NULL
    })
  }

  list(
    ncores = ncores,
    total_samples = total_samples,
    sample_ids = sample_ids,
    sample_chunks = sample_chunks,
    cl = cl
  )
}

#' @noRd
.h2_progress_update <- function(state, completed_units, total_units, label) {
  elapsed <- proc.time()[["elapsed"]] - state$start_time
  overall <- if (total_units > 0) completed_units / total_units else 1
  percent <- max(0, min(100, round(overall * 100)))
  msg <- paste0("[", sprintf("%3d", percent), "%] ", label)
  if (overall > 0 && overall < 1) {
    eta <- elapsed * (1 - overall) / overall
    eta <- max(0, round(eta))
    hh <- eta %/% 3600
    mm <- (eta %% 3600) %/% 60
    ss <- eta %% 60
    msg <- paste0(msg, " | Remaining ", sprintf("%02d:%02d:%02d", hh, mm, ss))
  } else if (overall >= 1) {
    elapsed <- max(0, round(elapsed))
    hh <- elapsed %/% 3600
    mm <- (elapsed %% 3600) %/% 60
    ss <- elapsed %% 60
    msg <- paste0(msg, " | total ", sprintf("%02d:%02d:%02d", hh, mm, ss))
  }

  if (!identical(msg, state$last_label)) {
    message(msg)
    state$last_label <- msg
  }
  invisible(NULL)
}
