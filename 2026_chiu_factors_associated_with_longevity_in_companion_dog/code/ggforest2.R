ggforest2 <- function (model, data = NULL, main = "Hazard ratio",
                       cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7,
                       refLabel = "reference", noDigits = 2,
                       xrange = NULL,
                       reference_free = FALSE,      # NEW: TRUE = deviations from overall mean
                       emm_weights = "cells")       # NEW: "cells" (counts) or "equal"
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  
  # --- coefficient table: standard (reference) OR reference-free (emmeans) ---
  if (!reference_free) {
    # ORIGINAL PATH: coefficients read directly off the model via tidy() [1]
    coef <- as.data.frame(tidy(model, conf.int = TRUE))
    rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
    
  } else {
    ## VERIFY #1 -----------------------------------------------------------
    ## reference_free assumes a SINGLE term of interest. If your model has
    ## more than one non-strata covariate, only the first is reparameterized.
    if (length(terms) != 1L)
      warning("reference_free assumes a single term of interest; using the first: ",
              names(terms)[1])
    rf_var <- names(terms)[1]
    
    em <- emmeans::emmeans(model, specs = rf_var, weights = emm_weights)
    
    ## CHANGED: effect contrasts = deviations from the GRAND MEAN across all
    ## levels. This gives EVERY level (incl. the former reference, e.g. WA)
    ## a genuine estimate, SE, CI, and p-value -- no level is pinned to 0.
    ## (Plain emmeans on a coxph returns log-HRs relative to the reference,
    ##  which zeroes WA; effect contrasts fix that.)
    emc   <- emmeans::contrast(em, method = "eff", infer = c(TRUE, TRUE))
    em.df <- as.data.frame(emc)
    
    ## VERIFY #2 -----------------------------------------------------------
    ## contrast() labels each row in a `contrast` column like "WA effect".
    ## Strip the " effect" suffix to recover the bare level ("WA") so that
    ## paste0(rf_var, lev) reproduces "stateWA" to match `inds` downstream.
    lev <- sub(" effect$", "", as.character(em.df$contrast))
    
    ## VERIFY #3 -----------------------------------------------------------
    ## contrast() on a coxph returns asymptotic CI names (asymp.LCL/UCL).
    ## Normalize to lower.CL/upper.CL so the data.frame() below is stable
    ## across emmeans variants.
    if (!"lower.CL" %in% names(em.df) && "asymp.LCL" %in% names(em.df)) {
      names(em.df)[names(em.df) == "asymp.LCL"] <- "lower.CL"
      names(em.df)[names(em.df) == "asymp.UCL"] <- "upper.CL"
    }
    
    coef <- data.frame(
      term      = paste0(rf_var, lev),   # must match `inds` format, see VERIFY #2
      estimate  = em.df$estimate,        # CHANGED: contrast() uses `estimate`, not `emmean`
      conf.low  = em.df$lower.CL,
      conf.high = em.df$upper.CL,
      ## VERIFY #4 -------------------------------------------------------
      ## p.value tests each level against the OVERALL MEAN (the centering
      ## point / dashed line at HR = 1). Stars mean "differs from the overall
      ## mean", NOT "differs from a reference". Reflect this in your caption.
      p.value   = if ("p.value" %in% names(em.df)) em.df$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
    rownames(coef) <- coef$term
    
    ## VERIFY #5 -----------------------------------------------------------
    ## No NA coefficients: effect contrasts return ALL levels (incl. WA) with
    ## estimates + CIs. The refLabel/"reference" logic below never triggers.
  }
  
  gmodel <- glance(model)
  
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (var %in% colnames(data)) {
      if (terms[i] %in% c("factor", "character")) {
        adf <- as.data.frame(table(data[, var]))   # N per level still correct [1]
        cbind(var = var, adf, pos = 1:nrow(adf))
      }
      else if (terms[i] == "numeric") {
        data.frame(var = var, Var1 = "", Freq = nrow(data), pos = 1)
      }
      else {
        vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
        data.frame(var = vars, Var1 = "", Freq = nrow(data), pos = seq_along(vars))
      }
    } else {
      message(var, "is not found in data columns, and will be skipped.")
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")   # e.g. "stateWA" [1]
  
  ## VERIFY #6 -------------------------------------------------------------
  ## Order alignment. allTermsDF is built from table(data[,var]), which is in
  ## FACTOR-LEVEL order; coef[inds, ] then reorders coef to match. As long as
  ## VERIFY #2 holds (level strings match), this join is correct even if
  ## emmeans returned rows in a different order. If toShow shows NAs, the
  ## mismatch is almost always the paste0 format in VERIFY #2.
  
  # ----- everything below is UNCHANGED from your original function [1] -----
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 0.05, 
                                                                  "*", ""), ifelse(toShowExpClean$p.value < 0.01, "*", 
                                                                                   ""), ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  if (is.null(xrange) | (length(xrange) > 2)) {
    rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                    na.rm = TRUE)
  } else if (length(xrange)==2) {
    rangeb <- log(xrange)
  } else {
    rangeb <- c(-1,1)*log(abs(xrange))
  }
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, 
                                                       "pt"), "mm"))

  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) + 
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]), 
                  fill = ordered(seq_along(var)%%2 + 1))) + scale_fill_manual(values = c("#FFFFFF33", 
                                                                                         "#00000033"), guide = "none") + geom_point(pch = 15, 
                                                                                                                                    size = 4) + geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), 
                                                                                                                                                              width = 0.15) + geom_hline(yintercept = 1, linetype = 3) + 
    coord_flip(ylim = exp(rangeplot)) + ggtitle(main) + scale_y_log10(name = "", 
                                                                      labels = sprintf("%g", breaks), expand = c(0.02, 0.02), 
                                                                      breaks = breaks) + theme_light() + theme(panel.grid.minor.y = element_blank(), 
                                                                                                               panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
                                                                                                               legend.position = "none", panel.border = element_blank(), 
                                                                                                               axis.title.y = element_blank(), axis.text.y = element_blank(), 
                                                                                                               axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
    xlab("") + annotate(geom = "text", x = x_annotate, y = exp(y_variable), 
                        label = toShowExpClean$var, fontface = "bold", hjust = 0, 
                        size = annot_size_mm) + annotate(geom = "text", x = x_annotate, 
                                                         y = exp(y_nlevel), hjust = 0, label = toShowExpClean$level, 
                                                         vjust = -0.1, size = annot_size_mm) + annotate(geom = "text", 
                                                                                                        x = x_annotate, y = exp(y_nlevel), label = toShowExpClean$N, 
                                                                                                        fontface = "italic", hjust = 0, vjust = ifelse(toShowExpClean$level == 
                                                                                                                                                         "", 0.5, 1.1), size = annot_size_mm) + annotate(geom = "text", 
                                                                                                                                                                                                         x = x_annotate, y = exp(y_cistring), label = toShowExpClean$estimate.1, 
                                                                                                                                                                                                         size = annot_size_mm, vjust = ifelse(toShowExpClean$estimate.1 == 
                                                                                                                                                                                                                                                "reference", 0.5, -0.1)) + annotate(geom = "text", 
                                                                                                                                                                                                                                                                                    x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci, 
                                                                                                                                                                                                                                                                                    size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    annotate(geom = "text", x = x_annotate, y = exp(y_stars), 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") + annotate(geom = "text", 
                                                           x = 0.5, y = exp(y_variable), label = paste0("# Events: ", 
                                                                                                        gmodel$nevent, "; Global p-value (Log-Rank): ", format.pval(gmodel$p.value.log, 
                                                                                                                                                                    eps = ".001"), " \nAIC: ", round(gmodel$AIC, 
                                                                                                                                                                                                     2), "; Concordance Index: ", round(gmodel$concordance, 
                                                                                                                                                                                                                                        2)), size = annot_size_mm, hjust = 0, vjust = 1.2, 
                                                           fontface = "italic")
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}
