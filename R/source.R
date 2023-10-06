.aprior <-
function (gamma.hat) 
{
    m <- mean(gamma.hat)
    s2 <- var(gamma.hat)
    (2 * s2 + m^2)/s2
}
.bprior <-
function (gamma.hat) 
{
    m <- mean(gamma.hat)
    s2 <- var(gamma.hat)
    (m * s2 + m^3)/s2
}
.combat_tmp1 <-
function (dat, batch, levels_batch, mod) 
{
    batchmod <- model.matrix(~-1 + batch)
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels_batch[i])
    }
    n.batches <- sapply(batches, length)
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    design <- as.matrix(design[, !check])
    batch.design <- design[, 1:n.batch]
    return(list(dat = dat, batchmod = batchmod, n.batch = n.batch, 
        batches = batches, n.batches = n.batches, n.array = n.array, 
        design = design, batch.design = batch.design))
}
.combat_tmp2 <-
function (tmp1, verbose = TRUE) 
{
    if (verbose) {
        cat("[combat] Adjusting for", ncol(tmp1$design) - ncol(tmp1$batchmod), 
            "covariate(s) or covariate level(s)\n")
    }
    if (qr(tmp1$design)$rank < ncol(tmp1$design)) {
        if (ncol(tmp1$design) == (tmp1$n.batch + 1)) {
            stop("[combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
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
    B.hat <- solve(t(tmp1$design) %*% tmp1$design) %*% t(tmp1$design) %*% 
        t(as.matrix(tmp1$dat))
    grand.mean <- t(tmp1$n.batches/tmp1$n.array) %*% B.hat[1:tmp1$n.batch, 
        ]
    var.pooled <- ((tmp1$dat - t(tmp1$design %*% B.hat))^2) %*% 
        rep(1/tmp1$n.array, tmp1$n.array)
    return(list(B.hat = B.hat, grand.mean = grand.mean, var.pooled = var.pooled))
}
.combat_tmp3 <-
function (dat, tmp1, tmp2, verbose = TRUE) 
{
    if (verbose) {
        cat("[combat] Standardizing data across features\n")
    }
    stand.mean <- t(tmp2$grand.mean) %*% t(rep(1, tmp1$n.array))
    if (!is.null(tmp1$design)) {
        tmp <- tmp1$design
        tmp[, c(1:tmp1$n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% tmp2$B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(tmp2$var.pooled) %*% t(rep(1, 
        tmp1$n.array)))
    return(list(stand.mean = stand.mean, s.data = s.data))
}
.combat_tmp4 <-
function (tmp1, tmp2, tmp3, eb = TRUE, verbose = TRUE) 
{
    if (eb) {
        if (verbose) {
            cat("[combat] Fitting L/S model and finding priors\n")
        }
    }
    else {
        if (verbose) {
            cat("[combat] Fitting L/S model\n")
        }
    }
    gamma.hat <- solve(t(tmp1$batch.design) %*% tmp1$batch.design) %*% 
        t(tmp1$batch.design) %*% t(as.matrix(tmp3$s.data))
    delta.hat <- NULL
    for (i in tmp1$batches) {
        delta.hat <- rbind(delta.hat, apply(tmp3$s.data[, i], 
            1, var, na.rm = T))
    }
    gamma.star <- delta.star <- NULL
    gamma.bar <- t2 <- a.prior <- b.prior <- NULL
    if (eb) {
        gamma.bar <- apply(gamma.hat, 1, mean)
        t2 <- apply(gamma.hat, 1, var)
        a.prior <- apply(delta.hat, 1, .aprior)
        b.prior <- apply(delta.hat, 1, .bprior)
        if (verbose) {
            cat("[combat] Finding parametric adjustments\n")
        }
        for (i in 1:tmp1$n.batch) {
            temp <- .it.sol(tmp3$s.data[, tmp1$batches[[i]]], 
                gamma.hat[i, ], delta.hat[i, ], gamma.bar[i], 
                t2[i], a.prior[i], b.prior[i])
            gamma.star <- rbind(gamma.star, temp[1, ])
            delta.star <- rbind(delta.star, temp[2, ])
        }
    }
    return(list(gamma.hat = gamma.hat, delta.hat = delta.hat, 
        gamma.star = gamma.star, delta.star = delta.star, gamma.bar = gamma.bar, 
        t2 = t2, a.prior = a.prior, b.prior = b.prior))
}
.combat_tmp5 <-
function (tmp1, tmp2, tmp3, tmp4, eb = TRUE, verbose = TRUE) 
{
    if (verbose) {
        cat("[combat] Adjusting the data\n")
    }
    bayesdata <- tmp3$s.data
    j <- 1
    for (i in tmp1$batches) {
        if (eb) {
            bayesdata[, i] <- (bayesdata[, i] - t(tmp1$batch.design[i, 
                ] %*% tmp4$gamma.star))/(sqrt(tmp4$delta.star[j, 
                ]) %*% t(rep(1, tmp1$n.batches[j])))
        }
        else {
            bayesdata[, i] <- (bayesdata[, i] - t(tmp1$batch.design[i, 
                ] %*% tmp4$gamma.hat))/(sqrt(tmp4$delta.hat[j, 
                ]) %*% t(rep(1, tmp1$n.batches[j])))
        }
        j <- j + 1
    }
    return((bayesdata * (sqrt(tmp2$var.pooled) %*% t(rep(1, tmp1$n.array)))) + 
        tmp3$stand.mean)
}
.it.sol <-
function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 0.0001) 
{
    n <- apply(!is.na(sdat), 1, sum)
    g.old <- g.hat
    d.old <- d.hat
    change <- 1
    count <- 0
    while (change > conv) {
        g.new <- .postmean(g.hat, g.bar, n, d.old, t2)
        sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
            1, sum, na.rm = T)
        d.new <- .postvar(sum2, n, a, b)
        change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
        g.old <- g.new
        d.old <- d.new
        count <- count + 1
    }
    adjust <- rbind(g.new, d.new)
    rownames(adjust) <- c("g.star", "d.star")
    adjust
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
combat_apply <-
function (combat_parameters, dat, site, cov = NULL, verbose = TRUE) 
{
    if (!is.factor(site) || nlevels(site) != length(combat_parameters$levels_batch) || 
        any(levels(site) != combat_parameters$levels_batch)) {
        stop("site must be a factor with the same levels than when fitting combat")
    }
    if (is.data.frame(dat)) {
        dat <- as.matrix(dat)
    }
    else if (!is.matrix(dat)) {
        stop("dat must be a matrix")
    }
    if (combat_parameters$transpose) {
        if (nrow(dat) != length(site)) {
            stop("dat must have the same number of rows than the length of site")
        }
        dat <- t(dat)
    }
    else {
        if (ncol(dat) != length(site)) {
            stop("dat must have the same number of columns than the length of site")
        }
    }
    dat.combat <- dat
    dat <- dat[combat_parameters$not_constant, ]
    if (is.data.frame((cov))) {
        cov <- as.matrix(cov)
    }
    else if (!(is.matrix(cov) || is.null(cov))) {
        stop("cov must be a matrix or NULL")
    }
    tmp1 <- .combat_tmp1(dat, site, combat_parameters$levels_batch, 
        cov)
    tmp3 <- .combat_tmp3(dat, tmp1, combat_parameters$tmp2, verbose)
    tmp5 <- .combat_tmp5(tmp1, combat_parameters$tmp2, tmp3, 
        combat_parameters$tmp4, combat_parameters$eb, verbose)
    dat.combat[combat_parameters$not_constant, ] <- tmp5
    if (combat_parameters$transpose) {
        dat.combat <- t(dat.combat)
    }
    return(list(dat.combat = dat.combat, gamma.hat = combat_parameters$tmp4$gamma.hat, 
        delta.hat = combat_parameters$tmp4$delta.hat, gamma.star = combat_parameters$tmp4$gamma.star, 
        delta.star = combat_parameters$tmp4$delta.star, gamma.bar = combat_parameters$tmp4$gamma.bar, 
        t2 = combat_parameters$tmp4$t2, a.prior = combat_parameters$tmp4$a.prior, 
        b.prior = combat_parameters$tmp4$b.prior, site = tmp1$batch, 
        cov = tmp1$cov, stand.mean = tmp3$stand.mean, stand.sd = sqrt(combat_parameters$tmp2$var.pooled)[, 
            1]))
}
combat_fit <-
function (dat, site, cov = NULL, eb = TRUE, verbose = TRUE) 
{
    if (!is.factor(site)) {
        stop("site must be a factor")
    }
    if (is.data.frame(dat)) {
        dat <- as.matrix(dat)
    }
    else if (!is.matrix(dat)) {
        stop("dat must be a matrix")
    }
    if (ncol(dat) == length(site)) {
        transpose <- FALSE
        if (verbose) {
            cat("[combat] Subjects are COLUMNS\n")
        }
    }
    else if (nrow(dat) == length(site)) {
        transpose <- TRUE
        dat <- t(dat)
        if (verbose) {
            cat("[combat] Subjects are ROWS\n")
        }
    }
    else {
        stop("dat must have the same number of columns or rows than the length of site")
    }
    if (is.data.frame((cov))) {
        cov <- as.matrix(cov)
    }
    else if (!(is.matrix(cov) || is.null(cov))) {
        stop("cov must be a matrix or NULL")
    }
    if (any(is.na(dat))) {
        if (verbose) {
            cat("[combat] Imputing missing data (only for fit)\n")
        }
        for (batch_i in sort(unique(site))) {
            i <- which(site == batch_i)
            dat_i <- dat[, i]
            if (any(is.na(dat_i))) {
                for (j in 1:nrow(dat)) {
                  dat_ji <- dat_i[j, ]
                  is_na <- which(is.na(dat_ji))
                  if (length(is_na) > 0 && length(is_na) < length(i)) {
                    if (length(is_na) == 1) {
                      mod_i_is_na <- matrix(cov[i[is_na], ], 
                        nrow = 1)
                    }
                    else {
                      mod_i_is_na <- cov[i[is_na], ]
                    }
                    beta <- matrix(coef(lm(dat_ji ~ cov[i, ])))
                    beta[which(is.na(beta))] <- 0
                    dat[j, i[is_na]] <- cbind(1, mod_i_is_na) %*% 
                      beta
                  }
                  else {
                    dat[j, i[is_na]] <- mean(dat_ji, na.rm = TRUE)
                  }
                }
            }
        }
        for (batch_i in sort(unique(site))) {
            i <- which(site == batch_i)
            if (any(is.na(dat[, i]))) {
                for (j in 1:nrow(dat)) {
                  dat_j <- dat[j, ]
                  if (is.na(dat_j[i[1]])) {
                    if (!is.null(cov)) {
                      beta <- matrix(coef(lm(dat_j ~ cov)))
                      beta[which(is.na(beta))] <- 0
                      dat[j, i] <- cbind(1, cov[i, ]) %*% beta
                    }
                    else {
                      dat[j, i] <- mean(dat_j, na.rm = TRUE)
                    }
                  }
                }
            }
        }
    }
    not_constant <- which(apply(dat, 1, function(x) {
        var(x) > 0
    }))
    dat <- dat[not_constant, ]
    if (eb) {
        if (verbose) {
            cat("[combat] Performing ComBat with empirical Bayes\n")
        }
    }
    else {
        if (verbose) {
            cat("[combat] Performing ComBat without empirical Bayes (L/S model)\n")
        }
    }
    if (verbose) {
        cat("[combat] Found", nlevels(site), "batches (e.g., sites)\n")
    }
    levels_batch <- levels(site)
    tmp1 <- .combat_tmp1(dat, site, levels_batch, cov)
    tmp2 <- .combat_tmp2(tmp1, verbose)
    tmp3 <- .combat_tmp3(dat, tmp1, tmp2, verbose)
    tmp4 <- .combat_tmp4(tmp1, tmp2, tmp3, eb, verbose)
    return(list(levels_batch = levels_batch, transpose = transpose, 
        not_constant = not_constant, eb = eb, tmp2 = tmp2, tmp4 = tmp4))
}
