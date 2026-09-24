# bendavid_2020b_hierarchical_vectors.R
#
# Specificity (negative agreement) and sensitivity (positive agreement) data
# from Bendavid et al. (2020b) Table in "Additional Data and Response to
# Comments" section. Source: Bendavid et al., medRxiv v2, April 30 2020.
# https://www.medrxiv.org/content/10.1101/2020.04.14.20062463v2.full.pdf
#
# Convention (Gelman & Carpenter 2020, book Ch 19):
#   index 1 = the study of interest (Bendavid site itself)
#   indices 2–13 = external specificity calibration studies
#   index 1 = the study of interest for sensitivity
#   indices 2–3 = external sensitivity calibration studies
#
# Specificity: y_spec = kit agreement N, n_spec = gold standard N
# Study 1 is the pooled manufacturer data used in Bendavid 2020a
# (index 1 represents spec[1] in the Stan model — the test's specificity
#  at the study site, informed by all calibration studies hierarchically)
#
# Reading from the Bendavid 2020b table (column order: Gold standard N,
# Kit agreement N):
#   1  371 / 368   Manufacturer pre-COVID sera
#   2   30 /  30   Pre-COVID healthy adults New York
#   3   70 /  70   Lyme-negative + early-2020 hospital
#   4 1102 / 1102  Pre-COVID hospital admissions
#   5  300 /  300  Pre-COVID plasma donors, heparin
#   6  311 /  311  Healthy adults, Biotest company staff
#   7  500 /  500  Children + adults, hospital, heparin
#   8  200 /  198  Pregnant women pre-COVID
#   9   99 /   99  Pre-COVID, RF < 450 U/ml
#  10   31 /   29  Pre-COVID, RF > 600 U/ml
#  11  150 /  146  China CDC, PCR-negative (COVID era)
#  12  108 /  105  Pre-COVID plasma donors
#  13   52 /   50  COVID-era PCR-negative, some other virus
#
# Note: Bendavid 2020b totals are 3324 / 3308, which matches summing below.
#   sum(n_spec) = 371+30+70+1102+300+311+500+200+99+31+150+108+52 = 3324 ✓
#   sum(y_spec) = 368+30+70+1102+300+311+500+198+99+29+146+105+50 = 3308 ✓

n_spec <- c(371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52)
y_spec <- c(368, 30, 70, 1102, 300, 311, 500, 198, 99,  29, 146, 105, 50)

# Sensitivity: y_sens = kit agreement N (positive tests), n_sens = gold std N
#   14   85 /  78  Manufacturer confirmed COVID patients (IgM data)
#   15   37 /  27  PCR-confirmed + IgG/IgM by ELISA
#   16   35 /  25  PCR-confirmed, 6-10 days from onset
#
# Note: Bendavid 2020b total 157 / 130.
#   sum(n_sens) = 85 + 37 + 35 = 157 ✓
#   sum(y_sens) = 78 + 27 + 25 = 130 ✓

n_sens <- c(85, 37, 35)
y_sens <- c(78, 27, 25)
