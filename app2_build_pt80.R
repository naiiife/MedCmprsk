## Adds the PT80 mediator to amp_survival_nab_analysis.csv (columns nab_max,
## ref_ic80_geo, pt80). Safe to rerun: previously added columns are dropped first.
##
## nab_max is the maximum predicted VRC01 serum concentration over study days 1-7
## (daily popPK predictions from Gilbert et al. 2022, Fig. 5 data). pt80 is nab_max
## over the geometric-mean IC80 (least sensitive isolate, right-censored at 100) of
## the viruses acquired by placebo-arm cases in the same trial. The file's own
## pt80_panel (nab_max2 / panel_ic80) is left unchanged.
library(dplyr)

nab_dir <- 'D:/Data/AMP/Neutralization Titer Biomarker/fig5_extFig7/datasets'
nab_map <- bind_rows(
  read.csv(file.path(nab_dir, 'id_uid_703.csv')) |> mutate(protocol = 'HVTN 703/HPTN 081'),
  read.csv(file.path(nab_dir, 'id_uid_704.csv')) |> mutate(protocol = 'HVTN 704/HPTN 085'))

nab7 <- read.csv(file.path(nab_dir, 'dat_est_daily_caseCtl.csv')) |>
  # one PK series per participant, except uid 9089 (enrolled in both trials), whose
  # two series are stacked; each restarts at time 1
  mutate(series = cumsum(time == 1), .by = uid) |>
  filter(time <= 7) |>
  summarise(nab_max = max(Cc), .by = c(uid, series)) |>
  inner_join(select(nab_map, uid, protocol), by = 'uid', relationship = 'many-to-many') |>
  # for 9089 the higher series is the 30 mg/kg HVTN 703 record, the lower the
  # 10 mg/kg HVTN 704 record
  filter(n() == 1 | (protocol == 'HVTN 703/HPTN 081') == (nab_max == max(nab_max)),
         .by = uid) |>
  select(uid, protocol, nab_max)
stopifnot(!anyDuplicated(nab7[c('uid', 'protocol')]))

ref_ic80 <- read.csv(file.path(nab_dir, 'dat_ic50_ic80.csv')) |>
  filter(dose == 'Control', !is.na(gmt80ls)) |>
  mutate(protocol = sub('HVTN(\\d+)/HPTN(\\d+)', 'HVTN \\1/HPTN \\2', study)) |>
  summarise(ref_ic80_geo = exp(mean(log(gmt80ls))), n_ref = n(), .by = protocol)

amp_file <- read.csv('amp_survival_nab_analysis.csv') |>
  select(-any_of(c('nab_max', 'ref_ic80_geo', 'pt80'))) |>
  left_join(nab7, by = c('uid', 'protocol')) |>
  left_join(select(ref_ic80, protocol, ref_ic80_geo), by = 'protocol') |>
  mutate(pt80 = nab_max / ref_ic80_geo)
stopifnot(all(amp_file$in_popPK_cohort[!is.na(amp_file$nab_max)]))

# readr writes full double precision; write.csv rounds to 15 significant digits
readr::write_csv(amp_file, 'amp_survival_nab_analysis.csv')
