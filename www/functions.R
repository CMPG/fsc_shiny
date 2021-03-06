### get entries labels
# takes a vector of pop size and returns entries labels
# Alexandre Gouy - 09-2017
getEntriesLabels <- function(pop.sizes) {
  
  if(prod(pop.sizes) >= 1e6) {
    stop("Too many entries (> 1e6).")
  }
  
  li <- lapply(pop.sizes, seq, from = 0, by = 1)
  entr <- expand.grid(li)
  entr <- rev(entr)
  entries.lab <- apply(entr, 1, function(x) {
    paste0("(", paste(x, collapse = ","), ")")
  })
  
  return(entries.lab)
}

# GETBARPLOTWORSTFITSFS
# function to look at the worst fitted entries of the joint SFS
# A. Gouy

plotWorstEntries <- function(obs, exp, n.entries = 30) {
  
  obs.SFS <- obs$sfs
  
  matrix.sfs <- getEntriesLabels(obs$pop.sizes)
  
  # Read and compute observed SFS, discarding entry for fixed ancestral and fixed derived
  obs$sfs[c(1, length(obs$sfs))] <- 0 
  
  exp.sfs <- exp$sfs * sum(obs$sfs)
  
  # compute and plot the relative difference between the SFS
  diff.sfs <- (obs$sfs - exp.sfs) / obs$sfs
  # par(mfrow = c(2, 1))
  # plot(diff.sfs, xaxt = 'n', type = "l", col = "black", main = "")
  # axis(1, at=seq(1, length(diff.sfs), by = 100), labels=matrix.sfs[seq(1,length(diff.sfs), by=100)], las=2, cex.axis=0.7)
  # 
  # compute the likelihood based on the expected sfs, but do this for each entry  
  eval <- obs$sfs > 0 & exp$sfs > 0
  exp.log.lik <- log10(exp$sfs[eval]) * obs$sfs[eval]  
  explhoodnomon <- sum(exp.log.lik, na.rm = T)    
  
  rel.obs.sfs <- obs$sfs / sum(obs$sfs)
  obs.log.lik <- log10(rel.obs.sfs[eval]) * obs$sfs[eval]
  obslhoodnomon <- sum(obs.log.lik, na.rm = T)
  
  # compute the difference in likelihood between expected and obs for each entry
  diff.log.lik <- exp.log.lik - obs.log.lik
  
  # correlation between the difference in likelihood and the number of sites observed
  # plot(obs$sfs[eval], diff.log.lik)
  
  # Get 30 worst entries with bad SFS fit
  outl <- order(abs(diff.log.lik), decreasing = TRUE)[seq_len(n.entries)]
  cat(matrix.sfs[eval][outl], sep="\n")
  
  par(
    mfrow = c(3, 1),
    mar = c(10, 6, 3, 3),
    las = 2
  ) 
  
  # get the entries with bias
  entries <- matrix.sfs[eval][outl]
  
  # plot the diff likelihood for those entries, based on relative SFS and on the observed counts
  rel_expSFS <- exp.sfs[match(entries, matrix.sfs)] / sum(exp.sfs)
  rel_obsSFS <- obs$sfs[match(entries, matrix.sfs)] / sum(obs$sfs)
  
  exp.log.lik <- log10(rel_expSFS) * obs$sfs[match(entries, matrix.sfs)]
  obs.log.lik <- log10(rel_obsSFS) * obs$sfs[match(entries, matrix.sfs)]
  diff.log.lik <- exp.log.lik - obs.log.lik  
  
  barplot(
    diff.log.lik,
    names.arg = entries,
    ylab = "Likelihood difference", 
    main = "Observed and expected SFS likelihoods differences"
  )
  
  # plot the barplots comparing counts for the expected and observed SFS (note that the expected SFS was multiplied by the number of SNPs)
  selentries <- match(entries, matrix.sfs)
  res <- matrix(
    c(exp.sfs[selentries], obs$sfs[selentries]),
    ncol = length(entries),
    byrow = TRUE
  )
  barplot(
    res,
    beside = TRUE,
    legend.text = c("Expected", "Observed"),
    names.arg = entries,
    # ylim = c(0, max(res) * 1.5),
    angle = 90,
    main = "Obs. and exp. SFS for worst entries"
  )
  
  # relative bias in the entries
  relres <- apply(res, 2, function(col){col[1]/col[2]})  
  barplot(
    relres,
    names.arg = entries,
    # ylim = c(0, max(relres) * 1.15),
    angle = 90,
    main = "Relative fit (expected / observed)",
    ylab="Relative fit"
  )
  abline(h = 1, lty = 2)
  
}


# COMPUTELHOOD
# computes the loglikelihood given an obs.sfs and exp.sfs
# INPUT
#   obs.sfs : numeric vector with the observed SFS
#   exp.sfs : numeric vector with the expected SFS
# OUTPUT
#   log likelihood computed as SUM m_i log10(p_i), 
#   where m_i is the ith entry of obs.sfs and p_i is the ith entry of exp.sfs
# NOTE: for entries where m_i > 0 and p_i=0, we replace p_i by a small value (penalty)
computelhood <- function(obs.sfs, exp.sfs) {
  
  lhood <- 0
  
  # remove the first and last entries
  obs.sfs <- obs.sfs[-c(1, length(obs.sfs))]
  exp.sfs <- exp.sfs[-c(1, length(exp.sfs))]
  
  # Get the valid entries, i.e. entries where obs.SFS > 0
  eval <- which(obs.sfs > 0)
  
  # Calculate expected SFS with the penaltie for entries where obs.SFS > 0 and exp.SFS == 0
  if(sum(exp.sfs[eval] == 0) > 0) {
    # Settings (penalty for exo SFS entries with zero)
    penalty <- 1e-10
    minpenalty <- 1e-8
    
    penalty <- min(exp.sfs[exp.sfs > 0]) / 100
    if(penalty > minpenalty) {
      penalty <- minpenalty
    } 
    
    # Get the entries which are zero at the obs SFS to have the penalty
    tmp.exp <- exp.sfs[eval]  # note that the length of tmp.exp is length(eval)
    tmp.exp[tmp.exp == 0] <- penalty 
    exp.sfs[eval] <- tmp.exp
  }
  
  # check that the sum of exp.sfs is 1
  exp.sfs <- exp.sfs / sum(exp.sfs)
  
  # compute the likelihood
  if(sum(exp.sfs[eval] == 0) > 0) {
    stop("Entries with exp.sfs = 0")
  } else {
    lhood <- sum(obs.sfs[eval] * log10(exp.sfs[eval]))    
  }
  
  lhood
}


getSFS <- function(filename) {
  
  obs.sfs <- readLines(con = filename, n = -1)
  entries <- as.numeric((strsplit(obs.sfs[[3]], "\t"))[[1]])
  pop.sizes <- as.numeric((strsplit(obs.sfs[[2]], "\t"))[[1]])[-1]
  
  npop <- length(pop.sizes)
  pop.names <- as.character(1:npop)
  
  return(list(sfs=entries,
              pop.sizes=pop.sizes,
              pop.names=pop.names,
              numpops=npop))
  
}

# Plot 1D SFS
plot1DSFS <- function(obs.SFS, exp.SFS, pop.sizes, pop.names, ylims, scaled = FALSE) {
  
  npops <- length(pop.sizes)
  if(missing(pop.names)) pop.names <- as.character(1:npops)
  
  obs.SFS[c(1,length(obs.SFS))] <- 0
  exp.SFS <- exp.SFS*sum(obs.SFS)
  exp.SFS[c(1,length(obs.SFS))] <- 0
  
  sfs.size <- pop.sizes+1  
  dim(obs.SFS) <- sfs.size[c(npops:1)]
  obs.SFS=aperm(obs.SFS, c(npops:1))
  
  dim(exp.SFS)=sfs.size[c(npops:1)]
  exp.SFS=aperm(exp.SFS, c(npops:1))
  
  # get the 1D marginals for each pop
  margobs <- list()
  margexp <- list()
  for(i in 1:npops) {
    margobs[[i]] <- apply(obs.SFS,i,sum)
    margexp[[i]] <- apply(exp.SFS,i,sum)
  }
  
  par(mfrow = c(1, npops), oma = c(2, 2, 5, 2))
  
  if(!scaled) {
    for(i in 1:npops) {
      plot(0:pop.sizes[i], log10(margobs[[i]]), 
           type = "s", pch = 2, ylim = ylims, 
           lwd = 2, xlab = pop.names[i],
           ylab = "log10(SNPcounts)")
      lines(0:pop.sizes[i], log10(margexp[[i]]), 
            type="s", pch=3, col=3, lty=1, lwd=2)
      legend("topright", c("Observed","Expected"), 
             col=c(1,3), lty=c(1,1), lwd=2)    
    }
    
    title(main="Observed vs. expected SFS", outer=T)
  }
  
  if(scaled) {
    for(i in 1:npops) {
      scaledobs <- margobs[[i]][-1]/sum(margobs[[i]][-1])
      scaledexp <- margexp[[i]][-1]/sum(margexp[[i]][-1])
      
      suri <- 1/(1:length(scaledobs))
      scal <- suri/sum(suri)
      
      statio <- 1/length(scaledobs)
      
      yobs <- statio * scaledobs / scal
      yexp <- statio * scaledexp / scal
      
      plot(1:pop.sizes[i], yobs, 
           type="s", pch=2, ylim=c(0, max(c(yobs, yexp))), 
           lwd=2, xlab=pop.names[i],
           ylab="log10(SNPcounts)")
      lines(1:pop.sizes[i], yexp, 
            type="s", pch=3, col=3, lty=1, lwd=2)
      abline(h = statio, lty = 2, col = 2)
      legend("topright", c("Observed","Expected", "Stationary"), 
             col=c(1, 3, 2), lty=c(1,1,2), lwd=2)    
    }
    
    title(main="Scaled observed vs. expected SFS", outer=T)
  }
  
}





read.input <- function(input.file) {
  
  ext <- substr(
    input.file, (nchar(input.file) + 1) - 3, nchar(input.file)
  )
  ## get input file(s)
  if(ext == "par") {
    fl <- input.file
  } else if(ext == "zip") {
    fl <- unzip(input.file)
  } else {
    stop("zip or par file expected.")
  }
  
  if(length(fl) == 0) {
    stop("Input file not found.")
  }
  
  # returns NA if not found:
  par.file <- fl[grep(".par", fl)[1]]
  obs.file <- fl[grep("DSFS.obs", fl)[1]]
  exp.file <- fl[grep("DSFS.txt", fl)[1]]
  lho.file <- fl[grep(".bestlhoods", fl)[1]]
  
  jsfs.file <- fl[grep("joint", fl)]
  
  f.names <- list(par = par.file,
                  obs = obs.file,
                  exp = exp.file,
                  lhood = lho.file,
                  jsfs = jsfs.file)
  
  return(f.names)
  
}



plot2DSFS <- function(input.file, freq.text = TRUE, freq.cutoff = 0.001, coul = "Red") {
  
  ext <- substr(
    input.file, (nchar(input.file) + 1) - 3, nchar(input.file)
  )
  
  if(ext == "obs") {
    joint <- read.delim(input.file, sep = "\t", header=TRUE, as.is=TRUE, skip = 1)
    rownames(joint) <- joint[, 1]
    joint.sfs <- as.matrix(joint[, -1])
    joint.sfs[1,1] <- 0
    joint.sfs <- joint.sfs/sum(joint.sfs)
  } else {
    joint <- read.delim(input.file, sep = "\t", header=TRUE, as.is=TRUE, skip = 0)
    rownames(joint) <- joint[, 1]
    joint.sfs <- as.matrix(joint[, -1])
  }

  ncut <- dim(joint.sfs)[1]
  
  if(coul == "Red") color <- rgb(0.8, 0.1, 0, seq(0, 1, length.out = ncut*10))
  else color <- suppressWarnings(RColorBrewer::brewer.pal(100, coul))
  
  image(log10(joint.sfs), col = color, asp = 0, bty = "n", xaxt = "n", yaxt = "n",
        xlab = "Population 1", ylab = "Population 2")
  axis(1, at = 1:ncut/ncut, labels = 1:ncut)
  axis(2, at = 1:ncut/ncut, labels = 1:ncut)
  
  if(freq.text) {
    coord <- expand.grid(0:(ncut-1)/(ncut-1), 0:(ncut-1)/(ncut-1))
    txt <- signif(joint.sfs, 1)
    txt[txt < freq.cutoff] <- NA
    text(coord[, 1], coord[, 2], txt, cex = 0.5)
  }

}


curvedArrows <- function (x1, y1, x2, y2, code = 2, size = 1, width = 1.2/4/cin, 
                          open = TRUE, sh.adj = 0.1, sh.lwd = 1, sh.col = if (is.R()) par("fg") else 1, 
                          sh.lty = 1, h.col = sh.col, h.col.bo = sh.col, h.lwd = sh.lwd, 
                          h.lty = sh.lty, curved = FALSE) 
{
  cin <- size * par("cin")[2]
  width <- width * (1.2/4/cin)
  uin <- if (is.R()) 
    1/xyinch()
  else par("uin")
  x <- sqrt(seq(0, cin^2, length = floor(35 * cin) + 2))
  delta <- sqrt(h.lwd) * par("cin")[2] * 0.005
  x.arr <- c(-rev(x), -x)
  wx2 <- width * x^2
  y.arr <- c(-rev(wx2 + delta), wx2 + delta)
  deg.arr <- c(atan2(y.arr, x.arr), NA)
  r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)
  bx1 <- x1
  bx2 <- x2
  by1 <- y1
  by2 <- y2
  lx <- length(x1)
  r.seg <- rep(cin * sh.adj, lx)
  theta1 <- atan2((y1 - y2) * uin[2], (x1 - x2) * uin[1])
  th.seg1 <- theta1 + rep(atan2(0, -cin), lx)
  theta2 <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
  th.seg2 <- theta2 + rep(atan2(0, -cin), lx)
  x1d <- y1d <- x2d <- y2d <- 0
  if (code %in% c(1, 3)) {
    x2d <- r.seg * cos(th.seg2)/uin[1]
    y2d <- r.seg * sin(th.seg2)/uin[2]
  }
  if (code %in% c(2, 3)) {
    x1d <- r.seg * cos(th.seg1)/uin[1]
    y1d <- r.seg * sin(th.seg1)/uin[2]
  }
  if (is.logical(curved) && all(!curved) || is.numeric(curved) && 
      all(!curved)) {
    segments(x1 + x1d, y1 + y1d, x2 + x2d, y2 + y2d, lwd = sh.lwd, 
             col = sh.col, lty = sh.lty)
    phi <- atan2(y1 - y2, x1 - x2)
    r <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
    lc.x <- x2 + 2/3 * r * cos(phi)
    lc.y <- y2 + 2/3 * r * sin(phi)
  }
  else {
    if (is.numeric(curved)) {
      lambda <- curved
    }
    else {
      lambda <- as.logical(curved) * 0.5
    }
    lambda <- rep(lambda, length.out = length(x1))
    c.x1 <- x1 + x1d
    c.y1 <- y1 + y1d
    c.x2 <- x2 + x2d
    c.y2 <- y2 + y2d
    midx <- (x1 + x2)/2
    midy <- (y1 + y2)/2
    spx <- midx - lambda * 1/2 * (c.y2 - c.y1)
    spy <- midy + lambda * 1/2 * (c.x2 - c.x1)
    sh.col <- rep(sh.col, length = length(c.x1))
    sh.lty <- rep(sh.lty, length = length(c.x1))
    sh.lwd <- rep(sh.lwd, length = length(c.x1))
    lc.x <- lc.y <- numeric(length(c.x1))
    for (i in seq_len(length(c.x1))) {
      if (lambda[i] == 0) {
        segments(c.x1[i], c.y1[i], c.x2[i], c.y2[i], 
                 lwd = sh.lwd[i], col = sh.col[i], lty = sh.lty[i])
        phi <- atan2(y1[i] - y2[i], x1[i] - x2[i])
        r <- sqrt((x1[i] - x2[i])^2 + (y1[i] - y2[i])^2)
        lc.x[i] <- x2[i] + 2/3 * r * cos(phi)
        lc.y[i] <- y2[i] + 2/3 * r * sin(phi)
      }
      else {
        spl <- xspline(x = c(c.x1[i], spx[i], c.x2[i]), 
                       y = c(c.y1[i], spy[i], c.y2[i]), shape = 1, 
                       draw = FALSE)
        lines(spl, lwd = sh.lwd[i], col = sh.col[i], 
              lty = sh.lty[i])
        if (code %in% c(2, 3)) {
          x1[i] <- spl$x[3 * length(spl$x)/4]
          y1[i] <- spl$y[3 * length(spl$y)/4]
        }
        if (code %in% c(1, 3)) {
          x2[i] <- spl$x[length(spl$x)/4]
          y2[i] <- spl$y[length(spl$y)/4]
        }
        lc.x[i] <- spl$x[2/3 * length(spl$x)]
        lc.y[i] <- spl$y[2/3 * length(spl$y)]
      }
    }
  }
  if (code %in% c(2, 3)) {
    theta <- atan2((by2 - y1) * uin[2], (bx2 - x1) * uin[1])
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(bx2, Rep)
    p.y2 <- rep(by2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)
    if (open) 
      lines((p.x2 + r.arr * cos(ttheta)/uin[1]), (p.y2 + 
                                                    r.arr * sin(ttheta)/uin[2]), lwd = h.lwd, col = h.col.bo, 
            lty = h.lty)
    else polygon(p.x2 + r.arr * cos(ttheta)/uin[1], p.y2 + 
                   r.arr * sin(ttheta)/uin[2], col = h.col, lwd = h.lwd, 
                 border = h.col.bo, lty = h.lty)
  }
  if (code %in% c(1, 3)) {
    x1 <- bx1
    y1 <- by1
    tmp <- x1
    x1 <- x2
    x2 <- tmp
    tmp <- y1
    y1 <- y2
    y2 <- tmp
    theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
    lx <- length(x1)
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(x2, Rep)
    p.y2 <- rep(y2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)
    if (open) 
      lines((p.x2 + r.arr * cos(ttheta)/uin[1]), (p.y2 + 
                                                    r.arr * sin(ttheta)/uin[2]), lwd = h.lwd, col = h.col.bo, 
            lty = h.lty)
    else polygon(p.x2 + r.arr * cos(ttheta)/uin[1], p.y2 + 
                   r.arr * sin(ttheta)/uin[2], col = h.col, lwd = h.lwd, 
                 border = h.col.bo, lty = h.lty)
  }
  list(lab.x = lc.x, lab.y = lc.y)
}

plotParFile <- function(
  input.file, gentime = 1, backwards = TRUE, inMar = 4.5,
  titulo = "Demographic model",
  growth = TRUE,
  cexMig = 0.5,
  cexAdm = 0.7,
  cexAx = 1,
  cexSide = 0.6,
  scalePop = 2.5,
  minPop = 0.001,
  arrowSplit = TRUE,
  plotMig = TRUE
) {
  
  ##### parameters
  codeArrow <- ifelse(backwards, 2, 1)
  par(mar = c(inMar+3, inMar, inMar, inMar), xpd = NA)
  
  ##### READ PAR FILE ------------------------------------------------------------
  parFile <- scan(
    input.file, what = character(0), sep = "\n",
    strip.white = TRUE, quiet = TRUE
  )
  ids <- grep("//", parFile)
  
  # get population sizes and number of samples
  pop.sizes <- as.numeric(parFile[(ids[2] + 1):(ids[3] - 1)])
  n.samp <- length(pop.sizes)
  
  if(n.samp == 1) {
    stop("Multiple populations are required to plot the parameter file.")
  }
  
  # get sampling times and sizes
  xu <- sapply(parFile[(ids[3] + 1):(ids[4] - 1)], strsplit, " ")
  samp.size <- as.numeric(unname(sapply(xu, "[", 1)))
  samp.times <- as.numeric(unname(sapply(xu, "[", 2)))
  samp.times[is.na(samp.times)] <- 0
  
  # get growth rates
  gr.rates <- as.numeric(parFile[(ids[4] + 1):(ids[5] - 1)])
  
  # get migration matrices
  nmig <- as.numeric(parFile[ids[5] + 1])
  
  if(nmig) {
    migMat <- list()
    for(i in 1:nmig) {
      migm <- parFile[(ids[5 + i] + 1):(ids[6 + i] - 1)]
      migm <- do.call("rbind", strsplit(migm, " "))
      class(migm) <- "numeric"
      migMat[[i]] <- migm
    }
  } else {
    migMat <- NA
  }
  
  # historical events
  histe <- parFile[(ids[6 + nmig] + 1):(ids[7 + nmig] - 1)][-1]
  if(length(histe) > 0){
    ma <- do.call("rbind", sapply(histe, strsplit, " "))
    if(dim(ma)[2] > 7) ma <- ma[, 1:7]
    rownames(ma) <- 1:length(histe)
    colnames(ma) <- c("time", "source", "sink", "mig", "newsize", "newgr", "migmat")
    suppressWarnings({class(ma) <- "numeric"})
    ma[, "source"] <- ma[, "source"] + 1
    ma[, "sink"] <- ma[, "sink"] + 1
  } else {
    ma <- NA
  }
  
  inp <- list(n.samp = n.samp, pop.sizes = pop.sizes,
              samp.size = samp.size, samp.times = samp.times,
              gr.rates = gr.rates, events = ma, migMat = migMat)
  
  nodes <- data.frame(id = 1:inp$n.samp, x = 1:inp$n.samp, y = inp$samp.times,
                      size = inp$pop.size, size.anc = NA, t.anc = NA, parent = NA, gr = inp$gr.rates)
  
  ##### PREPARE TREE -------------------------------------------------------------
  if(!is.na(inp$events[1])) {
    for(i in 1:dim(inp$events)[1]) {
      # if(inp$events[i, "mig"] == 1) {
      newid <- max(nodes$id) + 1
      newx <- inp$events[i, "sink"]
      newy <- inp$events[i, "time"]
      newsize <- inp$events[i, "newsize"] * nodes[inp$events[i, "sink"], "size"]
      
      nodes[c(inp$events[i, "source"]), "parent"] <- newid
      
      if(is.na(nodes[c(inp$events[i, "sink"]), "parent"])) {
        nodes[c(inp$events[i, "sink"]), "parent"] <- newid
      }
      nodes <- rbind(nodes, c(newid, newx, newy, newsize, NA, NA, NA, inp$events[i, "newgr"]))
      # }
    }
  }
  
  
  # add parents of new nodes
  for(i in 1:dim(nodes)[1]) {
    if(is.na(nodes[i, "parent"])) {
      nodes[i, "parent"] <- nodes[nodes[i, "x"], "parent"]
    }
  }
  
  # correct ancestral times if multiple nodes on the same x
  for(i in unique(nodes$x)) {
    nodes[nodes$x == i,]
    nl <- dim(nodes[nodes$x == i,])[1]
    if(nl > 1) {
      nodes[nodes$x == i,][-nl, ]$parent <- nodes[nodes$x == i,][-1, ]$id
    }
  }
  
  nodes$t.anc <- nodes[nodes$parent, "y"]
  nodes$size.anc <- nodes[nodes$parent, "size"]
  
  nodes$size.anc <- nodes$size*exp(nodes$gr*(nodes$t.anc-nodes$y))
  if(!is.na(inp$events[1])) maxT <- max(c(nodes$y, inp$events[, "time"]), na.rm = TRUE)
  else maxT <- max(nodes$y, na.rm = TRUE)
  
  plot(nodes$x, nodes$y, bty="n",
       col = "white", cex.axis = cexAx,
       xaxt = "n", xlab = "Deme",
       yaxt = "n", ylab = ifelse(gentime > 1, "Time (years)", "Time (generations)"),
       xlim = c(0, inp$n.samp + 1),
       ylim = c(0, maxT), main = titulo)
  axis(1, at = 1:inp$n.samp, labels = 1:inp$n.samp)
  axis(1, col = "white", tcl = 0, labels = NA)
  
  tmseq <- seq(0, signif(maxT, 2), signif(maxT, 2)/10)
  axis(2, at = tmseq, labels = tmseq * gentime)
  
  apply(nodes, 1, function(x) {
    scale <- max(c(nodes$size, nodes$size.anc))
    s0 <- x["size"] / scale
    s1 <- x["size.anc"] / scale
    s0 <- minPop + sqrt(s0) / scalePop
    s1 <- minPop + sqrt(s1) / scalePop
    
    ss <- 10^(1:5) / scale
    ss <- minPop + sqrt(ss) / scalePop
    
    segments(0, -seq(1,1.4,0.1) * (maxT * 0.15), 0+ss, -seq(1,1.4,0.1) * (maxT * 0.15),
             col = "bisque3", lwd = 4)
    text(0, -seq(1,1.4,0.1) * (maxT * 0.15), 10^(1:5),
         col = "black", lwd = 4, pos = 2, cex = 0.5)
    
    
    polygon(
      c(x["x"] - s0, x["x"] + s0, x["x"] + s1, x["x"]-s1),
      c(x["y"], x["y"], x["t.anc"], x["t.anc"]),
      col = "bisque3", border = NA
    )
    if(x["id"] == max(nodes$id)) {
      polygon(
        c(x["x"] - s0, x["x"] + s0, x["x"] + s0, x["x"]-s0),
        c(x["y"], x["y"], x["y"]*1.1, x["y"]*1.1),
        col = "bisque3", border = NA
      )
    }
  })
  
  if(!is.na(inp$events[1])) {
    for(i in 1:dim(nodes)[1]){
      if(nodes[nodes$parent[i], ]$x != nodes$x[i]) {
        if(arrowSplit) {
          curvedArrows(nodes$x[i], nodes$t.anc[i], nodes[nodes$parent[i], ]$x, nodes$t.anc[i],
                       sh.col = "bisque4", sh.lwd = 2, code = codeArrow, size = 0.5, width = 1)
        } else {
          segments(nodes$x[i], nodes$t.anc[i],
                   nodes[nodes$parent[i], ]$x, nodes$t.anc[i],
                   col = "bisque3", lwd = 2)
        }
      }
    }
    
    # introgression
    admi <- inp$events[0 < inp$events[, "mig"] & inp$events[, "mig"] < 1, ]
    if(length(admi) > 0) {
      if(length(admi) == 7) {
        curvedArrows(admi["source"], admi["time"], admi["sink"], admi["time"],
                     sh.col = "dodgerblue3", sh.lwd = 2, code = codeArrow, size = 0.5, width = 1)
        text(admi["sink"], admi["time"] - 0.02 * maxT,
             admi["mig"], col = "dodgerblue3",
             cex = cexAdm)
        
      } else {
        curvedArrows(admi[,"source"], admi[,"time"], admi[,"sink"], admi[,"time"],
                     sh.col = "dodgerblue3", sh.lwd = 2, code = codeArrow, size = 0.5, width = 1)
        text(admi[,"sink"], admi[,"time"] - 0.02 * maxT,
             admi[,"mig"], col = "dodgerblue3",
             cex = cexAdm)
        
      }
    }
    # events times
    tev <- inp$events[, "time"]
    text(inp$n.samp + 1, tev, tev, cex = cexSide, pos = c(2,4),
         offset = 0.2)
    segments(inp$n.samp + 1, 0, inp$n.samp + 1, max(inp$events[, "time"]), lty = 2)
    
    # migration matrix times
    mev <- inp$events[, "migmat"]
    text(0, tev, mev, cex = cexSide, pos = c(2,4),
         offset = 0.2)
    segments(0, 0, 0, max(inp$events[, "time"]), lty = 2)
  }
  
  # expansion/contraction
  if(growth) {
    typ <- nodes$gr; typ[typ<0] <- "^"; typ[typ>0] <- "v"; typ[typ==0] <- "="
    text(nodes$x, nodes$y, typ, col = "black", cex = 1, pos = 3, font = 2)
  }
  
  # migration
  if(nmig > 0 & plotMig){
    for(mi in 1:nmig) {
      if(mi == 1) t1 <- 0
      if(mi > 1) t1 <- inp$events[, "time"][inp$events[, "migmat"] == mi-1][1]
      m1 <- inp$migMat[[mi]]
      colnames(m1) <- 1:ncol(m1)
      rownames(m1) <- 1:nrow(m1)
      newm1 <- as.data.frame(as.table(m1))
      newm1 <- newm1[newm1$Freq > 0, ]
      
      if(nrow(newm1) > 0) {
        xm1 <- as.numeric(as.character(newm1[, 1]))
        xm2 <- as.numeric(as.character(newm1[, 2]))
        curvedArrows(xm1, rep(t1, nrow(newm1)), 
                     xm2, rep(t1, nrow(newm1)),
                     curve = maxT*0.05, sh.col = "forestgreen", size = 0.5, width = 1,
                     sh.lwd = 1)
        text((xm1 + xm2)/2, t1, signif(newm1[, 3], 2), pos = c(1, 3), cex = cexMig,
             col = "forestgreen")
      }
    }
  }
}