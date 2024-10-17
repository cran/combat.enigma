.aprior <-
function (gamma.hat) 
{
    m = mean(gamma.hat)
    s2 = var(gamma.hat)
    (2 * s2 + m^2)/s2
}
.bprior <-
function (gamma.hat) 
{
    m = mean(gamma.hat)
    s2 = var(gamma.hat)
    (m * s2 + m^3)/s2
}
.check_dat_site_and_cov <-
function (dat, site, cov, impute_missing_cov, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Checking dat\n")
    }
    if (!(is.matrix(dat) && is.numeric(dat))) {
        stop("dat must be a numeric matrix")
    }
    if (verbose) {
        cat("[combat.enigma] Checking site\n")
    }
    if (!is.factor(site)) {
        stop("site must be a factor")
    }
    if (any(table(site) == 0)) {
        stop("sites cannot be empty")
    }
    if (verbose) {
        cat("[combat.enigma] Checking cov\n")
    }
    if (!(impute_missing_cov %in% c(TRUE, FALSE))) {
        stop("impute_missing_cov must be TRUE or FALSE")
    }
    if (!(is.matrix(cov) && is.numeric(cov)) && !is.null(cov)) {
        stop("cov must be a numeric matrix or NULL")
    }
    if (any(is.na(cov)) & !impute_missing_cov) {
        stop("missing values in cov, use impute_missing_cov")
    }
}
.combat_tmp1 <-
function (dat, batch, levels_batch, mod) 
{
    batchmod = model.matrix(~-1 + batch)
    n.batch = nlevels(batch)
    batches = list()
    for (i in 1:n.batch) {
        batches[[i]] = which(batch == levels_batch[i])
    }
    n.batches = sapply(batches, length)
    n.array = sum(n.batches)
    design = cbind(batchmod, mod)
    check = apply(design, 2, function(x) all(x == 1))
    design = as.matrix(design[, !check])
    batch.design = design[, 1:n.batch]
    return(list(dat = dat, batchmod = batchmod, n.batch = n.batch, 
        batches = batches, n.batches = n.batches, n.array = n.array, 
        design = design, batch.design = batch.design))
}
.combat_tmp2 <-
function (tmp1, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Adjusting for", ncol(tmp1$design) - 
            ncol(tmp1$batchmod), "covariate(s) or covariate level(s)\n")
    }
    if (qr(tmp1$design)$rank < ncol(tmp1$design)) {
        if (ncol(tmp1$design) == (tmp1$n.batch + 1)) {
            stop("[combat.enigma] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
        }
        if (ncol(tmp1$design) > (tmp1$n.batch + 1)) {
            if ((qr(tmp1$design[, -c(1:tmp1$n.batch)])$rank < 
                ncol(tmp1$design[, -c(1:tmp1$n.batch)]))) {
                stop("The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.")
            }
            else {
                stop("At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
            }
        }
    }
    B.hat = solve(t(tmp1$design) %*% tmp1$design) %*% t(tmp1$design) %*% 
        t(as.matrix(tmp1$dat))
    grand.mean = t(tmp1$n.batches/tmp1$n.array) %*% B.hat[1:tmp1$n.batch, 
        ]
    var.pooled = ((tmp1$dat - t(tmp1$design %*% B.hat))^2) %*% 
        rep(1/tmp1$n.array, tmp1$n.array)
    return(list(B.hat = B.hat, grand.mean = grand.mean, var.pooled = var.pooled))
}
.combat_tmp3 <-
function (dat, tmp1, tmp2, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Standardizing data across features\n")
    }
    stand.mean = t(tmp2$grand.mean) %*% t(rep(1, tmp1$n.array))
    if (!is.null(tmp1$design)) {
        tmp = tmp1$design
        tmp[, c(1:tmp1$n.batch)] = 0
        stand.mean = stand.mean + t(tmp %*% tmp2$B.hat)
    }
    s.data = (dat - stand.mean)/(sqrt(tmp2$var.pooled) %*% t(rep(1, 
        tmp1$n.array)))
    return(list(stand.mean = stand.mean, s.data = s.data))
}
.combat_tmp4 <-
function (tmp1, tmp2, tmp3, eb, verbose) 
{
    if (eb) {
        if (verbose) {
            cat("[combat.enigma] Fitting L/S model and finding priors\n")
        }
    }
    else {
        if (verbose) {
            cat("[combat.enigma] Fitting L/S model\n")
        }
    }
    gamma.hat = solve(t(tmp1$batch.design) %*% tmp1$batch.design) %*% 
        t(tmp1$batch.design) %*% t(as.matrix(tmp3$s.data))
    delta.hat = NULL
    for (i in tmp1$batches) {
        delta.hat = rbind(delta.hat, apply(tmp3$s.data[, i], 
            1, var, na.rm = T))
    }
    gamma.star = delta.star = NULL
    gamma.bar = t2 = a.prior = b.prior = NULL
    if (eb) {
        gamma.bar = apply(gamma.hat, 1, mean)
        t2 = apply(gamma.hat, 1, var)
        a.prior = apply(delta.hat, 1, .aprior)
        b.prior = apply(delta.hat, 1, .bprior)
        if (verbose) {
            cat("[combat.enigma] Finding parametric adjustments\n")
        }
        for (i in 1:tmp1$n.batch) {
            temp = .it.sol(tmp3$s.data[, tmp1$batches[[i]]], 
                gamma.hat[i, ], delta.hat[i, ], gamma.bar[i], 
                t2[i], a.prior[i], b.prior[i])
            gamma.star = rbind(gamma.star, temp[1, ])
            delta.star = rbind(delta.star, temp[2, ])
        }
    }
    return(list(gamma.hat = gamma.hat, delta.hat = delta.hat, 
        gamma.star = gamma.star, delta.star = delta.star, gamma.bar = gamma.bar, 
        t2 = t2, a.prior = a.prior, b.prior = b.prior))
}
.combat_tmp5 <-
function (tmp1, tmp2, tmp3, tmp4, eb, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Adjusting the data\n")
    }
    bayesdata = tmp3$s.data
    j = 1
    for (i in tmp1$batches) {
        if (eb) {
            bayesdata[, i] = (bayesdata[, i] - t(tmp1$batch.design[i, 
                ] %*% tmp4$gamma.star))/(sqrt(tmp4$delta.star[j, 
                ]) %*% t(rep(1, tmp1$n.batches[j])))
        }
        else {
            bayesdata[, i] = (bayesdata[, i] - t(tmp1$batch.design[i, 
                ] %*% tmp4$gamma.hat))/(sqrt(tmp4$delta.hat[j, 
                ]) %*% t(rep(1, tmp1$n.batches[j])))
        }
        j = j + 1
    }
    return((bayesdata * (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))) + 
        tmp3$stand.mean)
}
.exclude_constant_voxels <-
function (dat, not_constant, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Excluding constant or NA voxels or regions\n")
    }
    dat[not_constant, ]
}
.exclude_small_sites_locally <-
function (dat, site, n.min, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Excluding sites with n <", n.min, 
            " locally (only for fit)\n")
    }
    if (any(table(site) < n.min)) {
        stop("At least one site has less subjects than n.min")
    }
    for (i in 1:nrow(dat)) {
        n_i = table(site[which(!is.na(dat[i, ]))])
        small_sites_i = names(n_i)[which(n_i < n.min & n_i > 
            0)]
        if (length(small_sites_i) > 0) {
            dat[i, which(site %in% small_sites_i)] = NA
        }
    }
    dat
}
.find_not_constant_voxels <-
function (dat, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Finding not constant or NA voxels or regions\n")
    }
    which(apply(dat, 1, function(x) {
        var(x, na.rm = TRUE) > 0
    }))
}
.find_prescaling_factors <-
function (dat, site, cov, levels_batch, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Finding prescaling factors\n")
    }
    dat_sd = NULL
    for (batch_i in levels_batch) {
        i = which(site == batch_i)
        if (!is.null(cov)) {
            cov_i = cov[i, ]
        }
        dat_sd_i = matrix(apply(dat[, i], 1, function(x) {
            if (any(!is.na(x))) {
                if (!is.null(cov)) {
                  sd(lm(x ~ cov_i)$residuals)
                }
                else {
                  sd(x)
                }
            }
            else {
                NA
            }
        }), nrow = 1)
        rownames(dat_sd_i) = batch_i
        dat_sd = rbind(dat_sd, dat_sd_i)
    }
    common_dat_sd = dat_sd[, which(apply(dat_sd, 2, function(x) {
        all(!is.na(x))
    }))]
    prescaling_factors = 1/apply(common_dat_sd, 1, median)
    if (verbose) {
        print(prescaling_factors)
    }
    prescaling_factors
}
.find_whether_dat_needs_transpose <-
function (dat, site, verbose) 
{
    if (ncol(dat) == length(site)) {
        if (verbose) {
            cat("[combat.enigma] Note: Subjects are COLUMNS\n")
        }
        FALSE
    }
    else if (nrow(dat) == length(site)) {
        if (verbose) {
            cat("[combat.enigma] Note: Subjects are ROWS\n")
        }
        TRUE
    }
    else {
        stop("dat must have the same number of columns or rows than the length of site")
    }
}
.impute_missing_cov <-
function (cov, site, cov_means, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Imputing missing cov\n")
    }
    for (batch_i in sort(unique(site))) {
        i = which(site == batch_i)
        cov_means_i = cov_means[which(rownames(cov_means) == 
            batch_i), ]
        for (j in 1:ncol(cov)) {
            is_na = which(is.na(cov[i, j]))
            if (length(is_na) > 0) {
                if (!is.na(cov_means_i[j])) {
                  cov[i[is_na], j] = cov_means_i[j]
                }
                else {
                  cov[i[is_na], j] = median(cov_means[, j], na.rm = TRUE)
                }
            }
        }
    }
    cov
}
.it.sol <-
function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) 
{
    n = apply(!is.na(sdat), 1, sum)
    g.old = g.hat
    d.old = d.hat
    change = 1
    count = 0
    while (change > conv) {
        g.new = .postmean(g.hat, g.bar, n, d.old, t2)
        sum2 = apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
            1, sum, na.rm = T)
        d.new = .postvar(sum2, n, a, b)
        change = max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
        g.old = g.new
        d.old = d.new
        count = count + 1
    }
    adjust = rbind(g.new, d.new)
    rownames(adjust) = c("g.star", "d.star")
    adjust
}
.lmm_tmp1 <-
function (dat, batch, levels_batch, mod, lin.hyp) 
{
    v_n.batch = c()
    m_random = NULL
    v_sigma = c()
    if (!is.null(mod)) {
        l_collinear = list()
    }
    l_coefficients = list()
    if (!is.null(lin.hyp)) {
        l_lin.hyp = list()
    }
    for (i in 1:nrow(dat)) {
        tryCatch({
            which_not_na_dat_i = which(!is.na(dat[i, ]))
            dat_i = dat[i, which_not_na_dat_i]
            batch_i = batch[which_not_na_dat_i]
            n.batch_i = length(unique(batch_i))
            if (!is.null(mod)) {
                mod_i = mod[which_not_na_dat_i, ]
                collinear_i = findLinearCombos(cbind(1, mod_i))$remove
                if (!is.null(collinear_i)) {
                  collinear_i = collinear_i - 1
                }
                constant_i = which(apply(t(t(mod_i)), 2, var) == 
                  0)
                if (length(constant_i) > 0) {
                  collinear_i = unique(c(collinear_i, constant_i))
                }
                if (!is.null(collinear_i)) {
                  mod_i = mod_i[, -collinear_i]
                }
            }
            if (!is.null(mod)) {
                if (ncol(t(t(mod_i))) > 0) {
                  if (n.batch_i > 1) {
                    m_i = lme(dat_i ~ mod_i, random = ~1 | batch_i)
                  }
                  else {
                    m_i = lm(dat_i ~ mod_i)
                  }
                }
                else {
                  if (n.batch_i > 1) {
                    m_i = lme(dat_i ~ 1, random = ~1 | batch_i)
                  }
                  else {
                    m_i = lm(dat_i ~ 1)
                  }
                }
            }
            else {
                if (n.batch_i > 1) {
                  m_i = lme(dat_i ~ 1, random = ~1 | batch_i)
                }
                else {
                  m_i = lm(dat_i ~ 1)
                }
            }
            summary_m_i = summary(m_i)
            v_n.batch = c(v_n.batch, n.batch_i)
            if (n.batch_i == 1) {
                random_i = matrix(0)
                rownames(random_i) = batch_i[1]
            }
            else {
                random_i = m_i$coefficients$random$batch
            }
            random_i = t(t(random_i[match(levels_batch, rownames(random_i)), 
                ]))
            rownames(random_i) = levels_batch
            colnames(random_i) = rownames(dat)[i]
            m_random = cbind(m_random, random_i)
            v_sigma = c(v_sigma, summary_m_i$sigma)
            if (!is.null(mod)) {
                if (is.null(collinear_i)) {
                  l_collinear[[i]] = NA
                }
                else {
                  l_collinear[[i]] = collinear_i
                }
            }
            if (n.batch_i == 1) {
                l_coefficients[[i]] = summary_m_i$coefficients
            }
            else {
                l_coefficients[[i]] = summary_m_i$tTable
            }
            if (!is.null(lin.hyp)) {
                if (is.null(collinear_i)) {
                  l_lin.hyp[[i]] = linearHypothesis(m_i, lin.hyp)
                }
                else if (rankMatrix(t(t(lin.hyp[, -collinear_i]))) == 
                  nrow(lin.hyp)) {
                  l_lin.hyp[[i]] = linearHypothesis(m_i, lin.hyp[, 
                    -collinear_i])
                }
                else {
                  l_lin.hyp[[i]] = NA
                }
            }
        }, error = function(e) {
            cat("Error in .lmm_tmp1, i =", i, "\n")
            cat("call:\n")
            print(e$call)
            cat("message:", e$message, "\n")
            browser()
        })
    }
    Out = list(n.batch = v_n.batch, random = m_random, sigma = v_sigma)
    if (!is.null(mod)) {
        Out$collinear = l_collinear
    }
    Out$coefficients = l_coefficients
    if (!is.null(lin.hyp)) {
        Out$lin.hyp = l_lin.hyp
    }
    return(Out)
}
.postmean <-
function (g.hat, g.bar, n, d.star, t2) 
{
    (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}
.postvar <-
function (sum2, n, a, b) 
{
    (0.5 * sum2 + b)/(n/2 + a - 1)
}
.prescale <-
function (dat, site, levels_batch, prescaling_factors, verbose) 
{
    if (verbose) {
        cat("[combat.enigma] Prescaling\n")
    }
    for (batch_i in levels_batch) {
        i = which(site == batch_i)
        dat[, i] = prescaling_factors[which(names(prescaling_factors) == 
            batch_i)] * dat[, i]
    }
    dat
}
combat_apply <-
function (combat_parameters, dat, site, cov = NULL, verbose = TRUE) 
{
    if (!(verbose %in% c(TRUE, FALSE))) {
        stop("verbose must be TRUE or FALSE")
    }
    if (!inherits(combat_parameters, "combat.enigma")) {
        stop("combat_parameters must be a list of class combat.enigma")
    }
    if (is.data.frame(dat)) {
        if (verbose) {
            cat("[combat.enigma] Converting dat into a matrix\n")
        }
        dat = as.matrix(dat)
    }
    if (is.data.frame((cov))) {
        cov = as.matrix(cov)
    }
    if (verbose) {
        cat("[combat.enigma] Checking dat\n")
    }
    if (!(is.matrix(dat) && is.numeric(dat))) {
        stop("dat must be a numeric matrix")
    }
    if (verbose) {
        cat("[combat.enigma] Checking site\n")
    }
    if (!is.factor(site) || nlevels(site) != length(combat_parameters$levels_batch) || 
        any(levels(site) != combat_parameters$levels_batch)) {
        stop("site must be a factor with the same levels than when fitting combat")
    }
    if (combat_parameters$transpose) {
        if (nrow(dat) != length(site)) {
            stop("site must have the same length as the number of rows of dat")
        }
    }
    else {
        if (ncol(dat) != length(site)) {
            stop("site must have the same length as the number of columns of dat")
        }
    }
    if (verbose) {
        cat("[combat.enigma] Checking cov\n")
    }
    if (!(is.matrix(cov) && is.numeric(cov)) && !is.null(cov)) {
        stop("cov must be a numeric matrix or NULL")
    }
    if (any(is.na(cov)) && !combat_parameters$impute_missing_cov) {
        stop("missing values in cov, use impute_missing_cov in combat_fit")
    }
    if (combat_parameters$transpose) {
        if (verbose) {
            cat("[combat.enigma] Transposing dat\n")
        }
        dat = t(dat)
    }
    if (combat_parameters$impute_missing_cov & any(is.na(cov))) {
        cov = .impute_missing_cov(cov, site, combat_parameters$cov_means, 
            verbose)
    }
    dat.combat = dat
    if (length(combat_parameters$not_constant) > 0) {
        dat = .exclude_constant_voxels(dat, combat_parameters$not_constant, 
            verbose)
    }
    if (combat_parameters$prescaling) {
        dat = .prescale(dat, site, combat_parameters$levels_batch, 
            combat_parameters$prescaling_factors, verbose)
    }
    if (combat_parameters$apply_combat) {
        combat_tmp1 = .combat_tmp1(dat, site, combat_parameters$levels_batch, 
            cov)
        combat_tmp3 = .combat_tmp3(dat, combat_tmp1, combat_parameters$combat_tmp2, 
            verbose)
        dat = .combat_tmp5(combat_tmp1, combat_parameters$combat_tmp2, 
            combat_tmp3, combat_parameters$combat_tmp4, combat_parameters$eb, 
            verbose)
    }
    dat.combat[combat_parameters$not_constant, ] = dat
    if (combat_parameters$transpose) {
        if (verbose) {
            cat("[combat.enigma] Back-transposing dat.combat\n")
        }
        dat.combat = t(dat.combat)
    }
    if (verbose) {
        cat("[combat.enigma] Note: combat_apply finished\n")
    }
    Out = list(dat.combat = dat.combat)
    if (combat_parameters$apply_combat) {
        Out$gamma.hat = combat_parameters$combat_tmp4$gamma.hat
        Out$delta.hat = combat_parameters$combat_tmp4$delta.hat
        Out$gamma.star = combat_parameters$combat_tmp4$gamma.star
        Out$delta.star = combat_parameters$combat_tmp4$delta.star
        Out$gamma.bar = combat_parameters$combat_tmp4$gamma.bar
        Out$t2 = combat_parameters$combat_tmp4$t2
        Out$a.prior = combat_parameters$combat_tmp4$a.prior
        Out$b.prior = combat_parameters$combat_tmp4$b.prior
        Out$stand.mean = combat_tmp3$stand.mean
        Out$stand.sd = sqrt(combat_parameters$combat_tmp2$var.pooled)[, 
            1]
    }
    Out
}
combat_fit <-
function (dat, site, cov = NULL, n.min = 10, impute_missing_cov = FALSE, 
    prescaling = FALSE, impute_empty_sites = FALSE, eb = TRUE, 
    verbose = TRUE) 
{
    if (!(verbose %in% c(TRUE, FALSE))) {
        stop("verbose must be TRUE or FALSE")
    }
    if (is.data.frame(dat)) {
        if (verbose) {
            cat("[combat.enigma] Converting dat into a matrix\n")
        }
        dat = as.matrix(dat)
    }
    if (is.data.frame((cov))) {
        cov = as.matrix(cov)
    }
    .check_dat_site_and_cov(dat, site, cov, impute_missing_cov, 
        verbose)
    levels_batch = levels(site)
    if (!(prescaling %in% c(TRUE, FALSE))) {
        stop("prescaling must be TRUE or FALSE")
    }
    if (!(impute_empty_sites %in% c(TRUE, FALSE))) {
        stop("impute_empty_sites must be TRUE or FALSE")
    }
    if (!(eb %in% c(TRUE, FALSE))) {
        stop("eb must be TRUE or FALSE")
    }
    transpose = .find_whether_dat_needs_transpose(dat, site, 
        verbose)
    if (transpose) {
        cat("[combat.enigma] Transposing dat\n")
        dat = t(dat)
    }
    if (impute_missing_cov) {
        cov_means = apply(cov, 2, by, site, mean, na.rm = TRUE)
        if (any(is.na(cov))) {
            cov = .impute_missing_cov(cov, site, cov_means, verbose)
        }
    }
    if (any(apply(dat, 1, function(dat_i) {
        min(table(site[which(!is.na(dat_i))]))
    }) < n.min)) {
        dat = .exclude_small_sites_locally(dat, site, n.min, 
            verbose)
    }
    not_constant = .find_not_constant_voxels(dat, verbose)
    if (length(not_constant) < nrow(dat)) {
        dat = .exclude_constant_voxels(dat, not_constant, verbose)
    }
    if (prescaling) {
        prescaling_factors = .find_prescaling_factors(dat, site, 
            cov, levels_batch, verbose)
        dat = .prescale(dat, site, levels_batch, prescaling_factors, 
            verbose)
    }
    if (any(is.na(dat))) {
        if (verbose) {
            cat("[combat.enigma] Imputing missing dat (only for fit)\n")
        }
        for (batch_i in levels_batch) {
            i = which(site == batch_i)
            dat_i = dat[, i]
            if (any(is.na(dat_i))) {
                for (j in 1:nrow(dat)) {
                  dat_ji = dat_i[j, ]
                  is_na = which(is.na(dat_ji))
                  if (length(is_na) > 0 && length(is_na) < length(i)) {
                    if (!is.null(cov)) {
                      if (length(is_na) == 1) {
                        mod_i_is_na = matrix(cov[i[is_na], ], 
                          nrow = 1)
                      }
                      else {
                        mod_i_is_na = cov[i[is_na], ]
                      }
                      beta = matrix(coef(lm(dat_ji ~ cov[i, ])))
                      beta[which(is.na(beta))] = 0
                      dat[j, i[is_na]] = cbind(1, mod_i_is_na) %*% 
                        beta
                    }
                    else {
                      dat[j, i[is_na]] = mean(dat_ji, na.rm = TRUE)
                    }
                  }
                }
            }
        }
        for (batch_i in levels_batch) {
            i = which(site == batch_i)
            if (any(is.na(dat[, i]))) {
                if (!impute_empty_sites) {
                  stop("empty site, use impute_empty_sites or lmm method")
                }
                for (j in 1:nrow(dat)) {
                  dat_j = dat[j, ]
                  if (is.na(dat_j[i[1]])) {
                    if (!is.null(cov)) {
                      beta = matrix(coef(lm(dat_j ~ cov)))
                      beta[which(is.na(beta))] = 0
                      dat[j, i] = cbind(1, cov[i, ]) %*% beta
                    }
                    else {
                      dat[j, i] = mean(dat_j, na.rm = TRUE)
                    }
                  }
                }
            }
        }
    }
    if (verbose) {
        if (eb) {
            cat("[combat.enigma] Fitting ComBat with empirical Bayes\n")
        }
        else {
            cat("[combat.enigma] Fitting ComBat without empirical Bayes\n")
        }
    }
    combat_tmp1 = .combat_tmp1(dat, site, levels_batch, cov)
    combat_tmp2 = .combat_tmp2(combat_tmp1, verbose)
    combat_tmp3 = .combat_tmp3(dat, combat_tmp1, combat_tmp2, 
        verbose)
    combat_tmp4 = .combat_tmp4(combat_tmp1, combat_tmp2, combat_tmp3, 
        eb, verbose)
    if (verbose) {
        cat("[combat.enigma] Note: combat_fit finished\n")
    }
    combat_parameters = list(apply_combat = TRUE, combat_tmp2 = combat_tmp2, 
        combat_tmp4 = combat_tmp4, eb = eb, impute_empty_sites = impute_empty_sites, 
        impute_missing_cov = impute_missing_cov, levels_batch = levels_batch, 
        not_constant = not_constant, prescaling = prescaling, 
        transpose = transpose)
    if (impute_missing_cov) {
        combat_parameters$cov_means = cov_means
    }
    if (prescaling) {
        combat_parameters$prescaling_factors = prescaling_factors
    }
    class(combat_parameters) = "combat.enigma"
    combat_parameters
}
lmm_fit <-
function (dat, site, cov = NULL, n.min = 10, impute_missing_cov = FALSE, 
    prescaling = FALSE, lin.hyp = NULL, verbose = TRUE) 
{
    if (!(verbose %in% c(TRUE, FALSE))) {
        stop("verbose must be TRUE or FALSE")
    }
    if (is.data.frame(dat)) {
        if (verbose) {
            cat("[combat.enigma] Converting dat into a matrix\n")
        }
        dat = as.matrix(dat)
    }
    if (is.data.frame((cov))) {
        cov = as.matrix(cov)
    }
    .check_dat_site_and_cov(dat, site, cov, impute_missing_cov, 
        verbose)
    levels_batch = levels(site)
    if (!(prescaling %in% c(TRUE, FALSE))) {
        stop("prescaling must be TRUE or FALSE")
    }
    transpose = .find_whether_dat_needs_transpose(dat, site, 
        verbose)
    if (transpose) {
        cat("[combat.enigma] Transposing dat\n")
        dat = t(dat)
    }
    initial_nrow_dat = nrow(dat)
    if (impute_missing_cov) {
        cov_means = apply(cov, 2, by, site, mean, na.rm = TRUE)
        if (any(is.na(cov))) {
            cov = .impute_missing_cov(cov, site, cov_means, verbose)
        }
    }
    if (any(apply(dat, 1, function(dat_i) {
        min(table(site[which(!is.na(dat_i))]))
    }) < n.min)) {
        dat = .exclude_small_sites_locally(dat, site, n.min, 
            verbose)
    }
    not_constant = .find_not_constant_voxels(dat, verbose)
    if (length(not_constant) < nrow(dat)) {
        dat = .exclude_constant_voxels(dat, not_constant, verbose)
    }
    if (prescaling) {
        prescaling_factors = .find_prescaling_factors(dat, site, 
            cov, levels_batch, verbose)
        dat = .prescale(dat, site, levels_batch, prescaling_factors, 
            verbose)
    }
    if (verbose) {
        cat("[combat.enigma] Fitting linear mixed-effects models\n")
    }
    lmm_tmp1 = .lmm_tmp1(dat, site, levels_batch, cov, lin.hyp)
    if (verbose) {
        cat("[combat.enigma] Creating results\n")
    }
    n.coef = 1 + ncol(cov)
    empty_map = rep(NA, initial_nrow_dat)
    b_map = t_map = z_map = list()
    for (coef_to_extract in 1:n.coef) {
        lmm_b = lmm_t = lmm_z = rep(NA, length(not_constant))
        for (i in 1:length(not_constant)) {
            if (length(lmm_tmp1$collinear[[i]]) == 1 && is.na(lmm_tmp1$collinear[[i]])) {
                b = lmm_tmp1$coefficients[[i]][coef_to_extract, 
                  ]
                lmm_b[i] = b[1]
                lmm_t[i] = b[length(b) - 1]
                lmm_z[i] = sign(b[1]) * (-qnorm(b[length(b)]/2))
            }
        }
        b_map[[coef_to_extract]] = empty_map
        b_map[[coef_to_extract]][not_constant] = lmm_b
        t_map[[coef_to_extract]] = empty_map
        t_map[[coef_to_extract]][not_constant] = lmm_t
        z_map[[coef_to_extract]] = empty_map
        z_map[[coef_to_extract]][not_constant] = lmm_z
    }
    sigma_map = empty_map
    sigma_map[not_constant] = lmm_tmp1$sigma
    if (!is.null(lin.hyp)) {
        lin.hyp_chisq = lin.hyp_z = rep(NA, length(not_constant))
        for (i in 1:length(not_constant)) {
            if (inherits(lmm_tmp1$lin.hyp[[i]], "anova")) {
                lin.hyp_chisq[i] = lmm_tmp1$lin.hyp[[i]]$Chisq[2]
                lin.hyp_z[i] = -qnorm(lmm_tmp1$lin.hyp[[i]]$`Pr(>Chisq)`[2])
            }
            else {
                if (!is.logical(lmm_tmp1$lin.hyp[[i]])) 
                  browser()
            }
        }
        lin.hyp_chisq_map = empty_map
        lin.hyp_chisq_map[not_constant] = lin.hyp_chisq
        lin.hyp_z_map = empty_map
        lin.hyp_z_map[not_constant] = lin.hyp_z
    }
    if (verbose) {
        cat("[combat.enigma] Note: lmm_fit finished\n")
    }
    combat_parameters = list(b_map = b_map, impute_missing_cov = impute_missing_cov, 
        levels_batch = levels_batch, lmm_tmp1 = lmm_tmp1, not_constant = not_constant, 
        prescaling = prescaling, sigma_map = sigma_map, t_map = t_map, 
        transpose = transpose, z_map = z_map)
    if (impute_missing_cov) {
        combat_parameters$cov_means = cov_means
    }
    if (prescaling) {
        combat_parameters$prescaling_factors = prescaling_factors
    }
    if (!is.null(lin.hyp)) {
        combat_parameters$lin.hyp_chisq_map = lin.hyp_chisq_map
        combat_parameters$lin.hyp_z_map = lin.hyp_z_map
    }
    class(combat_parameters) = "combat.enigma"
    combat_parameters
}
prescale_apply <-
function (combat_parameters, dat, site, cov = NULL, verbose = TRUE) 
{
    if (!(verbose %in% c(TRUE, FALSE))) {
        stop("verbose must be TRUE or FALSE")
    }
    if (!inherits(combat_parameters, "combat.enigma")) {
        stop("combat_parameters must be a list of class combat.enigma")
    }
    if (is.data.frame(dat)) {
        if (verbose) {
            cat("[combat.enigma] Converting dat into a matrix\n")
        }
        dat = as.matrix(dat)
    }
    if (is.data.frame((cov))) {
        cov = as.matrix(cov)
    }
    if (verbose) {
        cat("[combat.enigma] Checking dat\n")
    }
    if (!(is.matrix(dat) && is.numeric(dat))) {
        stop("dat must be a numeric matrix")
    }
    if (verbose) {
        cat("[combat.enigma] Checking site\n")
    }
    if (!is.factor(site) || nlevels(site) != length(combat_parameters$levels_batch) || 
        any(levels(site) != combat_parameters$levels_batch)) {
        stop("site must be a factor with the same levels than when fitting combat")
    }
    if (combat_parameters$transpose) {
        if (nrow(dat) != length(site)) {
            stop("site must have the same length as the number of rows of dat")
        }
    }
    else {
        if (ncol(dat) != length(site)) {
            stop("site must have the same length as the number of columns of dat")
        }
    }
    if (verbose) {
        cat("[combat.enigma] Checking cov\n")
    }
    if (!(is.matrix(cov) && is.numeric(cov)) && !is.null(cov)) {
        stop("cov must be a numeric matrix or NULL")
    }
    if (any(is.na(cov)) && !combat_parameters$impute_missing_cov) {
        stop("missing values in cov, use impute_missing_cov in combat_fit")
    }
    if (combat_parameters$transpose) {
        if (verbose) {
            cat("[combat.enigma] Transposing dat\n")
        }
        dat = t(dat)
    }
    if (combat_parameters$impute_missing_cov & any(is.na(cov))) {
        cov = .impute_missing_cov(cov, site, combat_parameters$cov_means, 
            verbose)
    }
    dat.combat = dat
    if (length(combat_parameters$not_constant) > 0) {
        dat = .exclude_constant_voxels(dat, combat_parameters$not_constant, 
            verbose)
    }
    if (combat_parameters$prescaling) {
        dat = .prescale(dat, site, combat_parameters$levels_batch, 
            combat_parameters$prescaling_factors, verbose)
    }
    if (combat_parameters$apply_combat) {
        combat_tmp1 = .combat_tmp1(dat, site, combat_parameters$levels_batch, 
            cov)
        combat_tmp3 = .combat_tmp3(dat, combat_tmp1, combat_parameters$combat_tmp2, 
            verbose)
        dat = .combat_tmp5(combat_tmp1, combat_parameters$combat_tmp2, 
            combat_tmp3, combat_parameters$combat_tmp4, combat_parameters$eb, 
            verbose)
    }
    dat.combat[combat_parameters$not_constant, ] = dat
    if (combat_parameters$transpose) {
        if (verbose) {
            cat("[combat.enigma] Back-transposing dat.combat\n")
        }
        dat.combat = t(dat.combat)
    }
    if (verbose) {
        cat("[combat.enigma] Note: combat_apply finished\n")
    }
    Out = list(dat.combat = dat.combat)
    if (combat_parameters$apply_combat) {
        Out$gamma.hat = combat_parameters$combat_tmp4$gamma.hat
        Out$delta.hat = combat_parameters$combat_tmp4$delta.hat
        Out$gamma.star = combat_parameters$combat_tmp4$gamma.star
        Out$delta.star = combat_parameters$combat_tmp4$delta.star
        Out$gamma.bar = combat_parameters$combat_tmp4$gamma.bar
        Out$t2 = combat_parameters$combat_tmp4$t2
        Out$a.prior = combat_parameters$combat_tmp4$a.prior
        Out$b.prior = combat_parameters$combat_tmp4$b.prior
        Out$stand.mean = combat_tmp3$stand.mean
        Out$stand.sd = sqrt(combat_parameters$combat_tmp2$var.pooled)[, 
            1]
    }
    Out
}
prescale_fit <-
function (dat, site, cov = NULL, n.min = 10, impute_missing_cov = FALSE, 
    verbose = TRUE) 
{
    if (!(verbose %in% c(TRUE, FALSE))) {
        stop("verbose must be TRUE or FALSE")
    }
    if (is.data.frame(dat)) {
        if (verbose) {
            cat("[combat.enigma] Converting dat into a matrix\n")
        }
        dat = as.matrix(dat)
    }
    if (is.data.frame((cov))) {
        cov = as.matrix(cov)
    }
    .check_dat_site_and_cov(dat, site, cov, impute_missing_cov, 
        verbose)
    levels_batch = levels(site)
    transpose = .find_whether_dat_needs_transpose(dat, site, 
        verbose)
    if (transpose) {
        cat("[combat.enigma] Transposing dat\n")
        dat = t(dat)
    }
    if (impute_missing_cov) {
        cov_means = apply(cov, 2, by, site, mean, na.rm = TRUE)
        if (any(is.na(cov))) {
            cov = .impute_missing_cov(cov, site, cov_means, verbose)
        }
    }
    if (any(apply(dat, 1, function(dat_i) {
        min(table(site[which(!is.na(dat_i))]))
    }) < n.min)) {
        dat = .exclude_small_sites_locally(dat, site, n.min, 
            verbose)
    }
    not_constant = .find_not_constant_voxels(dat, verbose)
    if (length(not_constant) < nrow(dat)) {
        dat = .exclude_constant_voxels(dat, not_constant, verbose)
    }
    prescaling_factors = .find_prescaling_factors(dat, site, 
        cov, levels_batch, verbose)
    dat = .prescale(dat, site, levels_batch, prescaling_factors, 
        verbose)
    if (verbose) {
        cat("[combat.enigma] Note: prescale_fit finished\n")
    }
    combat_parameters = list(apply_combat = FALSE, impute_missing_cov = impute_missing_cov, 
        levels_batch = levels_batch, not_constant = not_constant, 
        prescaling = TRUE, prescaling_factors = prescaling_factors, 
        transpose = transpose)
    if (impute_missing_cov) {
        combat_parameters$cov_means = cov_means
    }
    class(combat_parameters) = "combat.enigma"
    combat_parameters
}
