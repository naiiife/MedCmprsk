## Bootstrap for the frailty sensitivity analysis of application 1.
##
## Nonparametric bootstrap over participants of synth_cov2008.csv. Each replicate redraws
## the cross-fitting split (shared across the frailty values s, as in the main analysis)
## and refits everything, including the no-frailty total effect TE(0) that sets the scale
## factor -TE(0) / TE(s); fit_eif_cf_f with s = 0 reproduces fit_eif_cf exactly. The
## replicates are saved to app1_sens_boot.rds and read by run_sens() in
## application_analysis.R, which forms intervals as estimate +/- 1.96 x the robust (MAD)
## bootstrap SD.
##
## Run from the project folder; takes about 3 hours with 8 workers (B = 100, 2 events).

data_dir <- 'D:/Data/COVID_Vaccine'
B <- 100
s_seq <- seq(-1, 1, length.out = 11)
n_workers <- 8
out_file <- file.path(data_dir, 'app1_sens_boot.rds')

# Function definitions (and the quadrature globals) from application_analysis.R, without
# running any of its analyses
load_app_functions <- function(file = file.path(data_dir, 'application_analysis.R')) {
  env <- globalenv()
  for (e in parse(file)) {
    if (!is.call(e)) next
    head <- as.character(e[[1]])
    if (head == 'library') {
      eval(e, env)
    } else if (head %in% c('<-', '=') && is.name(e[[2]])) {
      rhs <- e[[3]]
      lhs <- as.character(e[[2]])
      if ((is.call(rhs) && identical(rhs[[1]], as.name('function'))) ||
          lhs %in% c('ngrid', 'quad', 'b_list', 'w_list')) eval(e, env)
    }
  }
  invisible(NULL)
}

boot_one <- function(b, event, s_seq, seed0 = 7316) {
  set.seed(seed0 + 1000 * event + b)
  d <- synth[sample.int(nrow(synth), replace = TRUE), ]
  n <- nrow(d)
  ss <- seq_len(n) %in% sample.int(n, floor(n / 2))
  do.call(rbind, lapply(s_seq, function(s) {
    out <- tryCatch({
      set.seed(seed0 + b)
      f <- fit_eif_cf_f(d$trt, d$Time, d$status, d$M, d[, Xname], method = 'riesz',
                        event = event, s = s, ss = ss, aug_time = FALSE)
      e <- est_decompose(f)
      c(TE = e$RMSTL_TE, ME = e$RMSTL_ME, DE = e$RMSTL_DE, CE = e$RMSTL_CE)
    }, error = function(err) c(TE = NA, ME = NA, DE = NA, CE = NA))
    data.frame(event = event, b = b, s = s, t(out))
  }))
}

if (sys.nframe() == 0L || identical(Sys.getenv('APP1_BOOT_RUN'), '1')) {
  library(parallel)
  cl <- makeCluster(n_workers)
  clusterExport(cl, c('data_dir', 'load_app_functions', 'boot_one'))
  clusterEvalQ(cl, {
    load_app_functions()
    synth <- read.csv(file.path(data_dir, 'synth_cov2008.csv'))
    Xname <- c('Vacc', 'Sex', 'Age')
    NULL
  })
  # one file per (event, replicate): progress can be monitored and an interrupted run
  # resumes where it stopped
  rep_dir <- file.path(data_dir, 'app1_sens_boot')
  dir.create(rep_dir, showWarnings = FALSE)
  jobs <- expand.grid(b = seq_len(B), event = 1:2)
  jobs$file <- file.path(rep_dir, sprintf('e%d_b%03d.rds', jobs$event, jobs$b))
  todo <- jobs[!file.exists(jobs$file), ]
  t0 <- Sys.time()
  parLapplyLB(cl, seq_len(nrow(todo)), function(i, todo, s_seq)
    saveRDS(boot_one(todo$b[i], todo$event[i], s_seq), todo$file[i]),
    todo = todo, s_seq = s_seq)
  stopCluster(cl)
  boot <- do.call(rbind, lapply(jobs$file, readRDS))
  attr(boot, 'elapsed') <- Sys.time() - t0
  saveRDS(boot, out_file)
}
