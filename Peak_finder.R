library(data.table)
library(RcppRoll)
########################################
# Calculate the pause sites on the positive strand.
# Import the Wig file containing the count of reads where 5' end locate on the positive strand.
original_wig_raw_plus = fread("./original_data_processed_sort_Positive.wig")
original_wig_plus = original_wig_raw_plus[, .(V2,V3)]
setnames(original_wig_plus,c("coordinate","count"))

setkey(original_wig_plus, coordinate)

original_wig_plus_median = original_wig_plus[ , Roll_median51 := roll_median(count, n = 51, weights = NULL, by = 1L, fill = NA, partial = FALSE,
                                                                   align = c("center", "left", "right"), normalize = TRUE, na.rm = FALSE),]
original_wig_plus_score = original_wig_plus_median[, score := count / Roll_median51]
original_wig_plus_socre_C = original_wig_plus_score[, score_F := score]
original_wig_plus_F = original_wig_plus_socre_C[score %in% c("Inf", "NaN"), score_F := count]
original_wig_plus_Final = original_wig_plus_F[, Strand := "+"]

# Export the pause sites meet the following two critetia: 1. Count >= 10/million reads; 2. Score >= 20.
# Reads_number represets the total number of reads uniquely aligned to the genome NC_000913.2.
pause_sites_plus = original_wig_plus_Final[(count >= (Reads_number/1000000)*10) & (score_F >= 20)]
########################################
original_wig_raw_minus = fread("../0.2 Bam files/Reads count and gene coverage analysis/SZ-6D2-1p5_S4_L012_R1_001_m14_M30_uniq_Negative.wig")

# Calculate the pause sites on the negative strand.
# Import the Wig file containing the count of reads where 5' end locate on the positive strand.
original_wig_raw_minus = fread("./original_data_processed_sort_Negative.wig")
original_wig_minus = original_wig_raw_minus[, .(V2,V3)]
setnames(original_wig_minus,c("coordinate","count"))

setkey(original_wig_minus, coordinate)

original_wig_minus_median = original_wig_minus[ , Roll_median51 := roll_median(count, n = 51, weights = NULL, by = 1L, fill = NA, partial = FALSE,
                                                                               align = c("center", "left", "right"), normalize = TRUE, na.rm = FALSE),]
original_wig_minus_score = original_wig_minus_median[, score := count / Roll_median51]
original_wig_minus_socre_C = original_wig_minus_score[, score_F := score]
original_wig_minus_F = original_wig_minus_socre_C[score %in% c("Inf", "NaN"), score_F := count]
original_wig_minus_Final = original_wig_minus_F[, Strand := "-"]

# Export the pause sites meet the following two critetia: 1. Count >= 10/million reads; 2. Score >= 20.
# Reads_number represets the total number of reads uniquely aligned to the genome NC_000913.2.
pause_sites_minus = original_wig_minus_Final[(count >= (Reads_number/1000000)*10) & (score_F >= 20)]
########################################
pause_sites = rbind(pause_sites_plus, pause_sites_minus)
