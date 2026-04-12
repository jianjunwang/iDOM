check_suggested <- function(pkg, caller) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required for %s(). Please install it and try again.",
        pkg,
        caller
      ),
      call. = FALSE
    )
  }
}

utils::globalVariables(
  c(
    ".",
    ".data",
    "Dist",
    "DK_all.degree",
    "DK_inner.degree",
    "DK_inter.degree",
    "Degree",
    "Degree.x",
    "Degree.y",
    "KK_all.degree",
    "KK_inner.degree",
    "KK_inter.degree",
    "Permutation",
    "Relative_Abundance",
    "Sample",
    "Trans.name",
    "Transformation_Database",
    "cDOM.IDs",
    "group",
    "iDME.inter",
    "iDME.intra",
    "m1",
    "m2",
    "n1",
    "n2",
    "num.trans.involved.in",
    "peak",
    "peak.x",
    "peak.y",
    "sd1",
    "sd2",
    "spearman.cor",
    "spearman.cor.p",
    "value",
    "variable"
  )
)
