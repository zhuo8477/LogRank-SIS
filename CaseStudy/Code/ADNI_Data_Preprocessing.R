library(tidyr)
library(ADNIMERGE2)
library(dplyr)
library(data.table)

data("DXSUM",   package = "ADNIMERGE2")
data("APOERES", package = "ADNIMERGE2")

baseline_mci <- DXSUM %>%
  filter(VISCODE == "bl", DIAGNOSIS == "MCI") %>%
  pull(PTID) %>%
  unique()

survival_Y <- DXSUM %>%
  filter(PTID %in% baseline_mci, !is.na(DIAGNOSIS)) %>%
  group_by(PTID) %>%
  mutate(BaseDate = min(EXAMDATE[VISCODE == "bl"], na.rm = TRUE)) %>%
  summarise(
    L_i    = max(as.numeric(EXAMDATE[DIAGNOSIS == "MCI"]      - BaseDate) / 30.44, na.rm = TRUE),
    R_i    = ifelse(
               any(DIAGNOSIS == "Dementia"),
               min(as.numeric(EXAMDATE[DIAGNOSIS == "Dementia"] - BaseDate) / 30.44, na.rm = TRUE),
               Inf
             ),
    status = as.integer(any(DIAGNOSIS == "Dementia"))
  ) %>%
  ungroup() %>%
  filter(L_i >= 0, L_i < R_i)

cat("Outcome data constructed:", nrow(survival_Y), "subjects\n")
cat("Event (status = 1):", sum(survival_Y$status == 1), "\n")
cat("Right-censored (status = 0):", sum(survival_Y$status == 0), "\n")

dt <- fread(file.choose())

cols_to_drop <- c("FID", "PAT", "MAT", "SEX", "PHENOTYPE")
dt[, (cols_to_drop) := NULL]
setnames(dt, "IID", "PTID")

ptid_list <- dt$PTID
dt[, PTID := NULL]
genomat           <- as.matrix(dt)
genomat[genomat == 2] <- 1L
rownames(genomat) <- ptid_list

common_ids <- intersect(survival_Y$PTID, rownames(genomat))
Y_final    <- survival_Y %>% filter(PTID %in% common_ids) %>% arrange(PTID)
X_final    <- genomat[Y_final$PTID, ]

cat("Final sample size  n =", nrow(Y_final), "\n")
cat("Number of SNPs     p =", ncol(X_final), "\n")

stopifnot(
  nrow(Y_final) == nrow(X_final),
  all(Y_final$PTID == rownames(X_final)),
  max(X_final, na.rm = TRUE) <= 1L,
  min(X_final, na.rm = TRUE) >= 0L,
  !any(Y_final$L_i >= Y_final$R_i)
)

cat("All validation checks passed.\n")
cat("Global minor allele frequency:", round(mean(X_final, na.rm = TRUE), 4), "\n")

saveRDS(Y_final, "ADNI_Y_Final.rds")
saveRDS(X_final, "ADNI_X_Final.rds")
