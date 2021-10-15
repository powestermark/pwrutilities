# Misc. initializers ------------------------------------------------------

#' Loads biocLite() from bioconductor.org
#'
#' @return NULL Called for side effect
#' @export
#'
sourcebioc <- function () {
  source("https://bioconductor.org/biocLite.R")
}



mygetmousemart <- function() {
  mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                           dataset="mmusculus_gene_ensembl",
                           host="www.ensembl.org")
  return(mart)
}



installed_date <- function (lib_index = 1L, mtime_first = FALSE) {
  if (Sys.info()["sysname"] == "Darwin" & !mtime_first) {
    tibble(package = list.files(.libPaths()[lib_index], full.names = TRUE)) %>%
      mutate(creation_date =
               lubridate::mdy_hms(
                 system(
                   str_c(c(
                     "echo",
                     package,
                     "| xargs -n1 GetFileInfo -d"
                   ), collapse = " "),
                   intern = TRUE))) %>%
      mutate(across(package, basename)) %>%
      arrange(desc(creation_date))
  } else {
    tibble(package = list.files(.libPaths()[lib_index], full.names = TRUE)) %>%
      mutate(file.info(package)) %>%
      dplyr::select(package, mtime) %>%
      mutate(across(package, basename)) %>%
      arrange(desc(mtime))
  }
}



# String utilities --------------------------------------------------------

returngoodmgi <- function(mgis) {
  ## assumes a vector, returns the "best" mgi symbol
  good <- grep("[0-9]Rik$", mgis, invert=TRUE)
  if (length(good) > 0) {
    return(mgis[good[1]])
  } else {
    return(mgis[1])
  }
}

weedriken <- function(mgi) {
  riks <- grepl("\\w+Rik", mgi)
  if (all(riks))
    ## return last one if only rikens
    return(mgi[length(mgi)])
  else {
    ## return last non-riken
    nonriks <- mgi[!riks]
    return(nonriks[length(nonriks)])
  }
}


weedmgi <- function (mgi, reference_mgi="xxxxxx") {
  ## if looking to choose a symbol from alias, first check whether a
  ## symbol is equal to a reference symbol
  if (any(mgi == reference_mgi))
    return(mgi[mgi == reference_mgi])
  ## Function to weed out RIKEN and such mgi symbols, if possible, and to keep
  ## only one mgi symbol in the end
  ## Start by sorting alphabetically, will return first hit (below)
  mgi <- sort(mgi)
  ## find out problematic mgi symbols
  riks <- grepl("\\w+Rik$", mgi)
  gms <- grepl("^Gm\\d+$", mgi)
  locs <- grepl("^LOC\\d+$", mgi)
  ## handle case when only problematic symbols are present
  if (all(riks | gms | locs)) {
    ## prefer returning Gm symbols
    if (any(gms))
      return(mgi[which(gms)[1]])
    else
      return(mgi[1])
  }
  ## otherwise, let's eliminate problematic symbols
  mgi <- mgi[!(riks | gms | locs)]
  return(mgi[1])
}


weed_refseq <- function (refseq_ids) {
  if (any(str_detect(refseq_ids, "^NM")))
    refseq_ids <- refseq_ids[str_detect(refseq_ids, "^NM")]
  if (any(str_detect(refseq_ids, "^NR")))
    refseq_ids <- refseq_ids[str_detect(refseq_ids, "^NR")]
  if (any(str_detect(refseq_ids, "^XM")))
    refseq_ids <- refseq_ids[str_detect(refseq_ids, "^XM")]
  if (any(str_detect(refseq_ids, "^XR")))
    refseq_ids <- refseq_ids[str_detect(refseq_ids, "^XR")]
  ## choose only one
  refseq_ids <- refseq_ids[which.min(nchar(refseq_ids))]
  refseq_ids
}


strrev <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")




#' Converts parameters to a string of par__value pairs
#'
#' @param ...
#'
#' @return string of par__value pairs separated by ___
#' @export
#'
#' @examples
#' a <- 1
#' b <- 2
#' var_value_string(a, b)
var_value_string <- function (...) {
  dots <- substitute(list(...))[-1]
  varnames <- sapply(dots, deparse)
  values <- mget(varnames, inherits=TRUE)
  paste(varnames, values, sep="__", collapse="___")
}



# Misc. utilities ---------------------------------------------------------

#' Round to nearest integer multiple of n
#'
#' @param x A vector
#' @param n A number
#'
#' @return Rounded vector
#' @export
#'
#' @examples
#' roundn(runif(10)*10, 2)
roundn <- function (x, n=1) {
  round(x/n)*n
}


#' Print the fraction of TRUEs in a vector of logicals
#'
#' @param vec Vector of logicals
#'
#' @return Fraction of TRUE values in vec
#' @export
#'
#' @examples
#' set.seed(123)
#' myfrac(rnorm(10) > 0)
#' myfrac(TRUE)
#' myfrac(FALSE)
myfrac <- function (vec) {
  tally <- table(vec)
  if (length(tally) == 2L)
    return(tally["TRUE"]/sum(tally))
  if (names(tally) == "FALSE")
    return(stats::setNames(0, "TRUE"))
  stats::setNames(1, "TRUE")
}

# curve of fraction true
myfraccurve <- function(pvec, limits)
  sapply(limits, function (x) myfrac(pvec < x))

val.field <- function(df, field) {
  vals <- df[, field]
  return(vals[!is.na(vals)])
}



#' Make a vector of stretches of indices of contigs satisfying a predicate
#' vector
#'
#' Useful together with split(), see example
#'
#' @param vec A vector
#' @param predicate A TRUE/FALSE vector of same length as vec, defaults to
#'   !is.na(vec)
#'
#' @return vector of indices
#' @export
#'
#' @examples
#' tosplit <- c(TRUE, TRUE, TRUE, NA, NA, TRUE, TRUE, NA,
#' TRUE, TRUE, TRUE, TRUE)
#' (index <- make_index_vector(tosplit))
#' split(tosplit[!is.na(tosplit)], index)
#'
#' set.seed(1234)
#' (tosplit2 <- rnorm(10))
#' make_index_vector(tosplit2, tosplit2 > 0)
#'
make_index_vector <- function (vec, predicate=!is.na(vec)) {
  vrle <- rle(predicate)
  ## values will be TRUE or FALSE
  vec_nonna_lengths <- vrle$lengths[vrle$values]
  rep(seq(1, length(vec_nonna_lengths)), times=vec_nonna_lengths)
}




#' Converts character vector to sorted factor
#'
#' Sorting is carried out using "human" criteria with gtools::mixedsort.
#' Useful with e.g., ggplot or tidyr::spread
#'
#' @param charvec A character vector
#' @param ... Further arguments to gtools::mixedsort
#'
#' @return Sorted factor
#' @export
#'
#' @examples
#' df <- dplyr::data_frame(cells=paste0("cell", 1:15), vals=1:15)
#' df %>% tidyr::spread(cells, vals)
#' df %>% dplyr::mutate(cells=char2sorted_factor(cells)) %>%
#' tidyr::spread(cells, vals)
#'
char2sorted_factor <- function(charvec, ...) {
  stopifnot(is.character(charvec))
  factor(charvec, levels=gtools::mixedsort(unique(charvec), ...))
}



#' Shifts vector by cyclic permutation
#'
#' @param vec vector to cycle
#' @param n amount of cyclic permutation
#'
#' @return shifted vector
#' @export
#'
#' @references
#' https://stackoverflow.com/questions/30542128/
#' circular-shifting-arrays-in-r-by-distance-n
#' @examples
#' shifter(1:5, 2)
#' shifter(1:5, -2)
shifter <- function (vec, n=1L) {
  if (n == 0) vec else c(tail(vec, -n), head(vec, n))
}




# Matrix utilities --------------------------------------------------------

returnslice <- function(position, vec, posminus, posplus) {
  return(vec[(position-posminus):(position+posplus)])
}

#' Matrix columnwise slice centered on different positions for each column
#'
#' @param mat A matrix
#' @param positions Center positions, vector. One element for each column in mat
#' @param posminus Scalar number of upstream positions in addition to center
#' @param posplus Scalar number of downstream positions in addition to center
#'
#' @return Matrix of column slices
#' @export
#'
#' @examples
#' A <- matrix(1:10, ncol=2)
#' matrixslices(A, c(4, 2), 1, 1)
#'
matrixslices <- function (mat, positions, posminus, posplus) {
  toret <- mapply(returnslice, positions,
                  split(mat, col(mat)),
                  MoreArgs=list(posminus=posminus, posplus=posplus))
  colnames(toret) <- colnames(mat)
  return(toret)
}

#' Running average
#'
#' @param vec Vector to perform running average on
#' @param winlength Averaging window length
#' @param keep.na Should NAs at beginning and end be kept?
#' (They start at the end)
#'
#' @return Running average
#' @export
#'
#' @examples
#' vec <- 1:4
#' running_average(vec, 2)
#' running_average(vec, 3)
#' running_average(vec, 2, keep.na=FALSE)
running_average <- function (vec, winlength, keep.na=TRUE) {
  runav <- stats::filter(vec, rep(1/winlength, winlength),
                         method="convolution") %>% as.vector()
  if (keep.na)
    runav
  else
    runav[!is.na(runav)]
}


compute.matrix.ratios <- function (mat1, mat2) {
  toret <- mat1/mat2
  toret[is.infinite(toret) | is.nan(toret)] <- NA
  return(toret)
}


m.head <- function (mat, n=6) mat[1:min(nrow(mat), n), 1:min(ncol(mat), n)]





# Circular helpers --------------------------------------------------------

circ_mean_24 <- function (x) {
  (as.numeric(circular::mean.circular(
    circular::circular(x, units="hours")))) %% 24
}

circ_sd_24 <- function (x) {
  as.numeric(circular::sd(circular::circular(x, units="hours")))*12/pi
}

circ_summary_24 <- function (x) {
  sumvec <- circular:::summary.circular(circular::circular(x, units="hours"))
  change_logic <- names(sumvec) %in% c("Min.", "1st Qu.", "Median",
                                       "Mean", "3rd Qu.", "Max.")
  sumvec[change_logic] <- sumvec[change_logic] %% 24
  sumvec
}

circ_kappa_24 <- function (x) {
  circular::circular(x, units="hours") %>%
    circular::mle.vonmises(bias=TRUE) %>%
    unclass %>%
    magrittr::extract(c("kappa", "se.kappa")) %>%
    unlist
}

circ_mu_24 <- function (x) {
  circular::circular(x, units="hours") %>%
    circular::mle.vonmises(bias=TRUE) %>%
    unclass %>%
    magrittr::extract(c("mu", "se.mu")) %>%
    unlist
}

## simple circular permutation around index
mycycleperm <- function (vec, ind) {
  if (!is.integer(ind)) ind <- as.integer(ind)
  if (ind == 1L)
    ## same phase
    return(vec)
  else
    ## reorder to put time matching nascent bin first
    return(vec[c(ind:length(vec), 1:(ind - 1))])
}

mycircpairdist.ct <- function(x, y) {
  delta.deg <- (x - y) %% 24
  delta.deg[delta.deg > 12] <- 24 - delta.deg[delta.deg > 12]
  return(delta.deg)
}

## alternative unsigned circular diff
ct.diff <- function (ct1, ct2) {
  cts <- cbind(ct1, ct2)
  return(apply(cts, 1, function(x) range(circular::circular(x, units="hours"))))
}

#' Signed differences in phase (0-24 hrs scale)
#'
#' Computes phase difference between x and y, between -12 and 12.
#' A positive result means that x has earlier phase than y. See example.
#'
#' @param x Vector, phases between 0 and 24
#' @param y Vector, phases between 0 and 24
#'
#' @return Vector of signed phase differences
#' @export
#'
#' @examples
#' mysignedcircpairdist.ct(c(10, 22, 2), c(12, 2, 23))
mysignedcircpairdist.ct <- function(x, y) {
  ## positive sign == x is ahead
  delta.deg <- (y - x)
  ## x more than 12h ahead == less than 12h behind
  delta.deg[delta.deg > 12 & !is.na(delta.deg)] <-
    -(24 - delta.deg[delta.deg > 12 & !is.na(delta.deg)])
  ## x more than 12h behind == less than 12h ahead
  delta.deg[delta.deg <= -12 & !is.na(delta.deg)] <-
    (24 + delta.deg[delta.deg <= -12 & !is.na(delta.deg)])
  return(delta.deg)
}

my.watson <- function(x, y) {
  xx <- circular::circular(x, units="hours")
  yy <- circular::circular(y, units="hours")
  circular::watson.two.test(xx, yy)
}


#'Split data frame containing ZTs over the 24/0 border
#'
#'@param df Data frame to split
#'@param grp Grouping variable (quosure), split will be made for each group
#'@param x Starting points, secondary variable in polar coordinates, e.g. length
#'  or amplitude
#'@param y Primary polar variable, e.g. time of day, ZT
#'@param xend Corresponding x end points
#'@param yend y end points
#'
#'@return Split data frame
#'@export
#'
#' @examples
#' # From plot in our paper 10.1073/pnas.1613103114
#' PFK <- 10.626
#' PFKup <- PFK + 2.9822/2
#' PFKdown <- PFK - 2.9822/2
#' ciopt <- dplyr::data_frame(x=1.25,
#'                     y=c(PFKdown + 7, PFKdown + 12) %% 24,
#'                         xend=1.25,
#'                         yend=c(PFKup + 7, PFKup + 12) %% 24,
#'                         what=c("PDH", "GS"))
#' ciopt
#' split024intervals(ciopt)
split024intervals <- function (df, grp=rlang::quo(what),
                               x="x", y="y",
                               xend="xend", yend="yend") {
  df %>% dplyr::group_by(!! grp) %>%
    dplyr::do (
      if (.[[yend]] > .[[y]])
        .
      else {
        tochange <- .[c(x, xend, y, yend)]
        therest <- setdiff(., tochange)
        changed <- setNames(dplyr::data_frame(x=.[[x]], xend=.[[xend]],
                                              y=c(.[[y]], 0),
                                              yend=c(24, .[[yend]])),
                            c(x, xend, y, yend))
        dplyr::bind_cols(dplyr::bind_rows(therest, therest), changed)
      }
    )
}


#' Constructs data_frame suitable for polar circadian phase arrows
#'
#' @param df input data_frame with:
#' @param amp_mean amplitudes (quosure or values overriding those in data frame)
#' @param phase_mean phases (quosure, see above)
#' @param grp grouping variable (quosure)
#'
#' @return data_frame suitable for ggplot
#' @export
#'
#' @examples
#' df <- dplyr::data_frame(what="test", amp=0.5, phase=12)
#' df %>% forge_polar_amp_df(rlang::quo(amp), rlang::quo(phase),
#' rlang::quo(what))
forge_polar_amp_df <- function (df, amp_mean, phase_mean,
                                grp=rlang::quo(what)) {
  df %>%
    dplyr::mutate(y = ((!! phase_mean)) %% 24,
                  yend = ((!! phase_mean)) %% 24,
                  x = 0, xend=(!! amp_mean)) %>%
    dplyr::select(!! grp, x, xend, y, yend)
}

#' Constructs data_frame suitable for polar circadian phase intervals
#'
#' @param df input data_frame with:
#' @param rad radial distance(s) for interval segment (quosure)
#' @param phase_mean phase means (quosure)
#' @param phase_int e.g. phase SD, positive (quosure)
#' @param grp grouping variable (quosure)
#'
#' @return data_frame suitable for ggplot
#' @export
#'
#' @examples
#' df <- dplyr::data_frame(what="test", rad=1, phase_mean=23, phase_int=2)
#' df %>% forge_polar_phaseint_df(rlang::quo(rad),
#' rlang::quo(phase_mean), rlang::quo(phase_int))
forge_polar_phaseint_df <- function (df, rad, phase_mean,
                                     phase_int,
                                     grp=rlang::quo(what)) {
  df %>%
    dplyr::mutate(y = ((!! phase_mean) - (!! phase_int)) %% 24,
                  yend = ((!! phase_mean) + (!! phase_int)) %% 24,
                  x = (!! rad), xend=(!! rad)) %>%
    dplyr::select(!! grp, x, xend, y, yend) %>%
    split024intervals(grp=grp)
}



# parallel cores ----------------------------------------------------------

mycores <- function() {
  return(min(30, parallel::detectCores()))
}



# ggplot: Basic publication themes ----------------------------------------

pw_theme <- ggplot2::theme_classic(base_size=7) +
  ggplot2::theme(
    #panel.background=element_rect(linetype="solid", color="black", size=1,
    #                              fill=NA),
    #panel.border=element_rect(linetype="solid",
    #                          color="black", size=1, fill="red"),
    #panel.grid.major.x=element_line(color="grey80"),
    #panel.spacing.y=unit(-0.5, "mm"),
    plot.margin=ggplot2::margin(3, 3, 3, 3),
    legend.key.size=ggplot2::unit(2.5, "mm"),
    legend.box.spacing=ggplot2::unit(0, "mm"),
    legend.justification=c(1, 1),
    legend.text=ggplot2::element_text(size=ggplot2::rel(1), color="black"),
    legend.title=ggplot2::element_text(size=ggplot2::rel(1), color="black"),
    axis.text=ggplot2::element_text(size=ggplot2::rel(1), color="black"),
    axis.line=ggplot2::element_line(size=0.5, lineend="round"),
    axis.ticks=ggplot2::element_line(color="black", lineend = "round")
  )

pw_theme_lattice <- ggplot2::theme(
  panel.background=ggplot2::element_rect(linetype="blank",
                                         fill=NA),
  panel.border=ggplot2::element_rect(linetype="solid",
                                     color="black", size=1, fill=NA),
  axis.line=ggplot2::element_blank(),
  strip.background=ggplot2::element_rect(fill="grey91", size=1),
  strip.placement="outside",
  strip.text=ggplot2::element_text(size=ggplot2::rel(1)),
  strip.text.x=ggplot2::element_text(margin=
                                       ggplot2::margin(2, 0, 2, 0, unit="pt")),
  strip.text.y=ggplot2::element_text(margin=
                                       ggplot2::margin(2, 0, 2, 0, unit="pt")),
  # panel.grid.major.x=element_line(color="grey80"),
  panel.spacing.y=ggplot2::unit(-0.5, "mm"),
)

pw_phase_theme <- pw_theme +
  ggplot2::theme(
    plot.margin=ggplot2::margin(0, 3, 0, 0, unit="pt"),
    panel.grid.major=ggplot2::element_line(color="grey92"),
    panel.grid.major.x=ggplot2::element_blank(),
    panel.grid.minor=ggplot2::element_blank(),
    panel.border=ggplot2::element_blank(),
    panel.spacing=ggplot2::unit(1, "cm"),
    plot.title=ggplot2::element_text(size=ggplot2::rel(1)),
    legend.margin=ggplot2::margin(0, 0, 0, 0, unit="lines"),
    legend.spacing=ggplot2::unit(0.6, "lines"),
    legend.key.width=ggplot2::unit(1, "lines"),
    legend.key.height=ggplot2::unit(0.6, "lines"),
    axis.line=ggplot2::element_blank(),
    axis.ticks=ggplot2::element_blank(),
    axis.text.y=ggplot2::element_blank(),
    axis.text.x=ggplot2::element_text(size=ggplot2::rel(1),
                                      color="black",
                                      margin=ggplot2::
                                        margin(0, 1, 0, 1, unit="cm")))

pw_theme_lattice_modern <-
  pw_theme_lattice +
  ggplot2::theme(
    legend.position = "top", legend.justification = "right",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "gray93", color = NA),
    panel.spacing.y = ggplot2::unit(1, "mm")
  )





# ggplot: Log 10 scales ---------------------------------------------------

## https://stackoverflow.com/questions/30179442/
## plotting-minor-breaks-on-a-log-scale-with-ggplot
log10_minor_break = function (...) {
  function (x) {
    minx         = floor(min(log10(x), na.rm=T)) - 1;
    maxx         = ceiling(max(log10(x), na.rm=T)) + 1;
    n_major      = maxx - minx + 1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks =
      rep(log10(seq(1, 9, by=1)), times=n_major) +
      rep(major_breaks, each=9)
    10^(minor_breaks)
  }
}

my_no0trail <- scales::format_format(drop0trailing=TRUE,
                                     scientific=FALSE)


my_x <- function (thebreaks=ggplot2::waiver(), ...) {
  ggplot2::scale_x_continuous(...,
                              breaks=thebreaks,
                              labels=my_no0trail)
}
my_y <- function (thebreaks=ggplot2::waiver(), ...) {
  ggplot2::scale_y_continuous(...,
                              breaks=thebreaks,
                              labels=my_no0trail)
}


my_xlog10 <- function (nbreaks=7,
                       thebreaks=scales::log_breaks(n=nbreaks), ...) {
  ggplot2::scale_x_log10(...,
                         breaks=thebreaks,
                         labels=scales::format_format(drop0trailing=TRUE,
                                                      scientific=FALSE))
}
my_ylog10 <- function (nbreaks=7,
                       thebreaks=scales::log_breaks(n=nbreaks), ...) {
  ggplot2::scale_y_log10(...,
                         breaks=thebreaks,
                         labels=scales::format_format(drop0trailing=TRUE,
                                                      scientific=FALSE))
}


my_xlog10p <- function (...) {
  ggplot2::scale_x_log10(breaks=scales::trans_breaks("log10",
                                                     function (x) 10^x,
                                                     ...),
                         labels=scales::
                           trans_format("log10", scales::math_format(10^.x)))
}
my_ylog10p <- function (...) {
  ggplot2::scale_y_log10(breaks=scales::trans_breaks("log10",
                                                     function (x) 10^x,
                                                     ...),
                         labels=scales::
                           trans_format("log10", scales::math_format(10^.x)))
}


## for e.g. kb
my_xlog10k <- function (nbreaks=7, scalefac=1000,
                        thebreaks=scales::log_breaks(n=nbreaks), ...) {
  ggplot2::scale_x_log10(breaks=thebreaks,
                         labels=scales::
                           trans_format(function (x) x/scalefac,
                                        scales::
                                          format_format(drop0trailing=TRUE,
                                                        scientific=FALSE)))
}
my_ylog10k <- function (nbreaks=7, scalefac=1000,
                        thebreaks=scales::log_breaks(n=nbreaks), ...) {
  ggplot2::scale_y_log10(breaks=thebreaks,
                         labels=scales::
                           trans_format(function (x) x/scalefac,
                                        scales::
                                          format_format(drop0trailing=TRUE,
                                                        scientific=FALSE)))
}



# ggplot: Circadian scales ------------------------------------------------

## returns function returning breaks with mystep intervals.
## limits gets rounded down and up, respectively, to multiples
## of mystep/rangedivfac.  There is a possibility of offset, to.
my_circbreaks <- function (mystep=12, rangedivfac=1, offset=0) {
  function (x) {
    myrangegrid <- mystep/rangedivfac
    minx <- x[1] - x[1] %% myrangegrid
    maxx <- x[2] + myrangegrid - x[2] %% myrangegrid
    seq(minx, maxx, by=mystep) + offset
  }
}

## usage example
# my_x_circscale <- function (mystep=12, rangedivfac=1, ...) {
#   thebreaks <- my_circbreaks(mystep, rangedivfac)
#   scale_x_continuous(breaks=thebreaks,
#                      limits=c(thebreaks[1], thebreaks[length(thebreaks)]),
#                      ...)
# }
#
# my_y_circscale <- function (mystep=12, rangedivfac=1, ...) {
#   thebreaks <- my_circbreaks(mystep, rangedivfac)
#   scale_y_continuous(breaks=thebreaks,
#                      limits=c(thebreaks[1], thebreaks[length(thebreaks)]),
#                      ...)
# }



# ggplot: circadian extensions --------------------------------------------


log_relamp_2_relamp <- function (means, log_relamp)
  (2^(2*means*log_relamp) - 1)/(2^(2*means*log_relamp) + 1)

log_amp_2_relamp <- function (log_amp)
  (2^(2*log_amp) - 1)/(2^(2*log_amp) + 1)


compute_hrfit <- function (data, scales, params, n=101L, per=24) {

  rng <- range(data$x, na.rm = TRUE)
  grid <- data.frame(x = seq(rng[1], rng[2], length=n))
  logdata <- data
  logdata$y <- log2(data$y)

  mod <- Rfit::rfit(y ~ cos(2*pi/per*x) + sin(2*pi/per*x),
                    data = logdata)

  rmean_l <- coef(mod)[1]
  rcos_l <- coef(mod)[2]
  rsin_l <- coef(mod)[3]
  amp_l <- sqrt(rcos_l^2 + rsin_l^2)
  phi_l <- atan2(rsin_l, rcos_l) %% (2*pi)

  rmean <- 2^rmean_l
  relamp <- log_amp_2_relamp(amp_l)

  grid$y <- rmean + rmean*relamp*cos(2*pi/per*grid$x - phi_l)
  grid
}



# _________________



compute_fit_dlnorm <- function (data, scales, params,
                                na.rm = TRUE, n = 101L,
                                freq = TRUE, binwidth = 1) {

  rng <- range(data$x, na.rm = na.rm)
  grid <- data.frame(x = seq(rng[1], rng[2], length = n))
  grid$y <- dlnorm(grid$x, meanlog = mean(log(data$x)), sdlog = sd(log(data$x)))

  if(freq) {
    grid$y <- grid$y*length(data$x)*binwidth
  }

  grid
}


StatFitDlnorm <- ggplot2::ggproto("StatFitDlnorm", ggplot2::Stat,
                                  required_aes = "x",
                                  compute_group = compute_fit_dlnorm)


stat_fit_dlnorm <- function (mapping = NULL, data = NULL, geom = "line",
                             position = "identity", na.rm = TRUE,
                             show.legend = NA, inherit.aes = TRUE,
                             n = 101L, freq = TRUE, binwidth = 1, ...) {
  ggplot2::layer(
    stat = StatFitDlnorm, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, freq = freq, binwidth = binwidth, ...)
  )
}


# _________________



compute_fit_dnorm <- function (data, scales, params,
                               na.rm = TRUE, n = 101L,
                               freq = TRUE, binwidth = 1) {

  rng <- range(data$x, na.rm = na.rm)
  grid <- data.frame(x = seq(rng[1], rng[2], length = n))
  grid$y <- dnorm(grid$x, mean = mean(data$x), sd = sd(data$x))

  if(freq) {
    grid$y <- grid$y*length(data$x)*binwidth
  }

  grid
}


StatFitDnorm <- ggplot2::ggproto("StatFitDnorm", ggplot2::Stat,
                                 required_aes = "x",
                                 compute_group = compute_fit_dnorm)


stat_fit_dnorm <- function (mapping = NULL, data = NULL, geom = "line",
                            position = "identity", na.rm = TRUE,
                            show.legend = NA, inherit.aes = TRUE,
                            n = 101L, freq = TRUE, binwidth = 1, ...) {
  ggplot2::layer(
    stat = StatFitDnorm, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, freq = freq, binwidth = binwidth, ...)
  )
}


# _________________



compute_pfit_logtrans <- function (data, scales, params,
                                   na.rm = TRUE, n = 101L) {

  rng <- range(data$x, na.rm = na.rm)
  grid <- data.frame(x = seq(rng[1], rng[2], length = n))
  logdata <- data
  logdata$x <- log(data$x)
  logdata$y <- log(data$y)

  # bootstrap weights for weighted linear least squares
  boot_fit <- lm(y ~ x, data = logdata)
  # obtain averaged standardized variance as function of y
  residual_fit <- loess(rstandard(boot_fit)^2 ~ fitted(boot_fit))
  # use inverse variance as weights for the final weighted fit
  fit <- lm(y ~ x, data = logdata, weights = 1/fitted(residual_fit))

  grid$y <- exp(predict(fit, data.frame(x = log(grid$x))))

  grid
}


StatPfitLog <- ggplot2::ggproto("StatPfitLog", ggplot2::Stat,
                                required_aes = c("x", "y"),
                                compute_group = compute_pfit_logtrans)


stat_pfit <- function (mapping = NULL, data = NULL, geom = "line",
                       position = "identity", na.rm = TRUE, show.legend = NA,
                       inherit.aes = TRUE, n = 101L, ...) {
  ggplot2::layer(
    stat = StatPfitLog, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, ...)
  )
}



# _________________



compute_hrfit_logtrans <- function (data, scales, params, na.rm = TRUE,
                                    n = 101L, per = 24) {

  rng <- range(data$x, na.rm = na.rm)
  grid <- data.frame(x = seq(rng[1], rng[2], length=n))
  logdata <- data
  logdata$y <- log2(data$y)

  mod <- Rfit::rfit(y ~ cos(2*pi/per*x) + sin(2*pi/per*x),
                    data = logdata)
  ## there's no predict() in Rfit
  grid$y <- 2^(coef(mod)[1] + coef(mod)[2]*cos(2*pi/per*grid$x) +
                 coef(mod)[3]*sin(2*pi/per*grid$x))
  grid
}


StatHrfitLog <- ggplot2::ggproto("StatHrfitLog", ggplot2::Stat,
                                 required_aes = c("x", "y"),
                                 compute_group = compute_hrfit_logtrans)


stat_hrfit <- function (mapping = NULL, data = NULL, geom = "line",
                        position = "identity", na.rm = TRUE, show.legend = NA,
                        inherit.aes = TRUE, n=101L, per=24, ...) {
  ggplot2::layer(
    stat = StatHrfitLog, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, n = n, per = per, ...)
  )
}



# _________________



compute_hrhline_logtrans <- function (data, scales, params,
                                      na.rm = TRUE, per = 24) {

  rng <- range(data$x, na.rm = na.rm)
  grid <- data.frame(x = c(rng[1], rng[2]))
  logdata <- data
  logdata$y <- log2(data$y)

  mod <- Rfit::rfit(y ~ cos(2*pi/per*x) + sin(2*pi/per*x),
                    data = logdata)
  ## there's no predict() in Rfit
  grid$y <- rep(2^(coef(mod)[1]), 2)
  grid
}


StatHrhlineLog <- ggplot2::ggproto("StatHrhlineLog", ggplot2::Stat,
                                   required_aes = c("x", "y"),
                                   compute_group = compute_hrhline_logtrans)


stat_hrhline <- function (mapping = NULL, data = NULL, geom = "line",
                          position = "identity", na.rm = TRUE,
                          show.legend = NA,
                          inherit.aes = TRUE, per=24, ...) {
  ggplot2::layer(
    stat = StatHrhlineLog, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, per = per, ...)
  )
}



# _________________



compute_density_phase <- function (data, scales, params, bw=6) {

  densest <-
    circular::density.circular(circular::circular(data$x, units="hour"), bw=bw)
  grid <- data.frame(x = unclass(densest$x),
                     y = densest$y*pi/12)

  grid
}

StatDensPhase <- ggplot2::ggproto("StatDensPhase", ggplot2::Stat,
                                  required_aes="x",
                                  compute_group=compute_density_phase)

stat_density_phase <- function (mapping = NULL, data = NULL, geom = "line",
                                position = "identity", na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE, bw = 6, ...) {
  ggplot2::layer(
    stat = StatDensPhase, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm=na.rm, bw=bw, ...)
  )
}


# ggplot: other extensions ------------------------------------------------


## uperrorbar from
## https://stackoverflow.com/questions/27585776/
## error-bars-for-barplot-only-in-one-direction

geom_uperrorbar <- function (mapping = NULL, data = NULL,
                             stat = "identity", position = "identity",
                             ...,
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


GeomUperrorbar <-
  ggplot2::ggproto("GeomUperrorbar", ggplot2::Geom,
                   default_aes = ggplot2::aes(colour = "black", size = 0.5,
                                              linetype = 1, width = 0.5,
                                              alpha = NA),
                   draw_key = ggplot2::draw_key_path,
                   required_aes = c("x", "y", "ymax"),
                   setup_data = function (data, params) {
                     data$width <- data$width %||%
                       params$width %||%
                       (ggplot2::resolution(data$x, FALSE) * 0.9)
                     transform(data,
                               xmin = x - width / 2,
                               xmax = x + width / 2, width = NULL
                     )
                   },
                   draw_panel = function (data, panel_scales, coord,
                                          width = NULL) {
                     ggplot2::GeomPath$draw_panel(data.frame(
                       x = as.vector(rbind(data$xmin, data$xmax, NA,
                                           data$x,   data$x)),
                       y = as.vector(rbind(data$ymax, data$ymax, NA,
                                           data$ymax, data$y)),
                       colour = rep(data$colour, each = 5),
                       alpha = rep(data$alpha, each = 5),
                       size = rep(data$size, each = 5),
                       linetype = rep(data$linetype, each = 5),
                       group = rep(1:(nrow(data)), each = 5),
                       stringsAsFactors = FALSE,
                       row.names = 1:(nrow(data) * 5)
                     ), panel_scales, coord)
                   }
  )


## NOTE
## can then use geom="uperrorbar" in stat_summary


# ggplot: custom plots ----------------------------------------------------

my_hist <- function (values, n = 10, breaks = NA,
                     binwidth = NA, boundary = floor(min(values)),
                     breaks_x = 6, breaks_y = 10,
                     color = "black", fill = "grey95") {
  if (!all(is.na(breaks))) {
    the_hist <- ggplot2::geom_histogram(
      breaks = breaks,
      color = color, fill = fill
    )
  } else if (!is.na(binwidth)) {
    the_hist <- ggplot2::geom_histogram(
      binwidth = binwidth, boundary = boundary,
      color = color, fill = fill
    )
  } else {
    the_hist <- ggplot2::geom_histogram(
      bins = n,
      color = color, fill = fill
    )
  }

  if (is.function(breaks_x))
    thebreaks_x <- breaks_x
  else
    thebreaks_x <- scales::breaks_extended(breaks_x)

  if (is.function(breaks_y))
    thebreaks_y <- breaks_y
  else
    thebreaks_y <- scales::breaks_extended(breaks_y)

  dplyr::tibble(x = values) %>%
    ggplot2::ggplot(ggplot2::aes(x)) +
    the_hist +
    my_x(expand = ggplot2::expand_scale(mult = 0.01),
         thebreaks = thebreaks_x) +
    my_y(expand = ggplot2::expand_scale(mult = c(0.002, 0.05)),
         thebreaks = thebreaks_y) +
    pw_theme +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90"),
      panel.grid.minor.y = ggplot2::element_line(color = "grey90")
    )
}


# Graphics ----------------------------------------------------------------


#' Computes color components along a circular spectrum
#'
#' Uses RColorBrewer as basis
#'
#' @param ncolor Number of desired colors
#' @param ncolor.orig Number of colors in the original spectrum, defaults to 11,
#' the maximal number of available colors in RColorBrewer
#'
#' @return Desired colors in hex format
#' @export
#'
#' @examples
#' graphics::pie(rep(1, 24), col=circular_spectral_colors(24L))
circular_spectral_colors <- function (ncolor=24L, ncolor.orig=11L) {
  all_cb_cols <- RColorBrewer::brewer.pal(ncolor.orig, "Spectral")
  ## get a cyclically defined ramp function
  ## It has to go all the way, otherwise a big gap between the endings will
  ## arise.
  cyclic_ramp_fun <- grDevices::colorRampPalette(c(all_cb_cols, all_cb_cols[1]))
  ## return requested amount of colors, let the ramp go all the way around then
  ## discard the last value
  cyclic_ramp_fun(ncolor + 1)[1:ncolor]
}



logexp2plotmath <- function(theseq) {
  sapply(theseq,
         function(x) as.expression(substitute(10^y, list(y=x))))
}




# Genomics ----------------------------------------------------------------

#' List of genes with introns
#'
#' @param txdb TxDb database object
#'
#' @return GRangesList with flattened genes with introns
#' @export
#'
#' @references See post on (https://stat.ethz.ch/pipermail/bioconductor/
#' 2013-February/050725.html)
#'
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm9.refGene)
#' my_intronsByGene(TxDb.Mmusculus.UCSC.mm9.refGene)
#'
my_intronsByGene <- function (txdb) {
  introns <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
  ## this creates a GRanges with all introns, not ordered according to anything
  ulst <- unlist(introns)
  ## for now, just remove dupes
  intronsNoDups <- ulst[!duplicated(ulst)]

  ## obtains a gene id for each transcript name
  suppressMessages(
    themap <- AnnotationDbi::select(txdb,
                                    keys=names(intronsNoDups),
                                    keytype="TXNAME", columns="GENEID"))
  ## there must be no NAs
  stopifnot(all(!is.na(themap$GENEID)))
  ## this works only when there is a 1:1 or many:1 mapping between txid and
  ## geneid
  if (length(intronsNoDups) != nrow(themap))
    stop(paste("Intron to gene annotation is not possible in this annotation,",
               "one-to-many mapping(s) between transcript and gene present"))

  S4Vectors::split(intronsNoDups, themap$GENEID)
}

#' List of genes with introns
#'
#' @param txdb TxDb database object
#'
#' @return GRangesList with flattened genes with introns
#' @export
#'
#' @references See post on (https://stat.ethz.ch/pipermail/
#' bioconductor/2013-February/050725.html)
#'
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm9.refGene)
#' my_intronsByGene(TxDb.Mmusculus.UCSC.mm9.refGene)
#'
my_intronsByGene2 <- function (txdb) {
  introns <- GenomicFeatures::intronsByTranscript(txdb, use.names=TRUE)
  ## this creates a GRanges with all introns, not ordered according to anything
  ulst <- unlist(introns)
  ## for now, just remove dupes
  intronsNoDups <- ulst[!duplicated(ulst)]

  classic_names <- names(intronsNoDups) %>%
    stringr::str_extract("^\\w+")
  names(intronsNoDups) <- classic_names

  ## obtains a gene id for each transcript name
  suppressMessages(
    themap <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                    keys=classic_names,
                                    keytype="REFSEQ", columns="ENTREZID"))
  ## there must be no NAs
  stopifnot(all(!is.na(themap$ENTREZID)))
  ## this works only when there is a 1:1 or many:1 mapping between txid and
  ## geneid
  if (length(intronsNoDups) != nrow(themap))
    stop(paste("Intron to gene annotation is not possible in this annotation,",
               "one-to-many mapping(s) between transcript and gene present"))

  S4Vectors::split(intronsNoDups, themap$ENTREZID)
}

#' List of genes with exons
#'
#' @param txdb TxDb database object
#'
#' @return GRangesList with flattened genes with introns
#' @export
#'
#' @references See post on (https://stat.ethz.ch/pipermail/
#' bioconductor/2013-February/050725.html)
#'
#' @examples
#' library(TxDb.Mmusculus.UCSC.mm9.refGene)
#' my_intronsByGene(TxDb.Mmusculus.UCSC.mm9.refGene)
#'
my_exonsByGene2 <- function (txdb) {
  exons <- GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)
  ## this creates a GRanges with all introns, not ordered according to anything
  ulst <- unlist(exons)
  ## for now, just remove dupes
  exonsNoDups <- ulst[!duplicated(ulst)]

  classic_names <- names(exonsNoDups) %>%
    stringr::str_extract("^\\w+")
  names(exonsNoDups) <- classic_names

  ## obtains a gene id for each transcript name
  suppressMessages(
    themap <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                    keys=classic_names,
                                    keytype="REFSEQ", columns="ENTREZID"))
  ## there must be no NAs
  stopifnot(all(!is.na(themap$ENTREZID)))
  ## this works only when there is a 1:1 or many:1 mapping between txid and
  ## geneid
  if (length(exonsNoDups) != nrow(themap))
    stop(paste("Exon to gene annotation is not possible in this annotation,",
               "one-to-many mapping(s) between transcript and gene present"))

  S4Vectors::split(exonsNoDups, themap$ENTREZID)
}


m2e <- function (mgisyms) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=mgisyms,
                              keytype="SYMBOL",
                              columns="ENTREZID")
  df$ENTREZID
}

e2m <- function (entrezids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=entrezids,
                              keytype="ENTREZID",
                              columns="SYMBOL")
  df$SYMBOL
}

e2me <- function (entrezids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=entrezids,
                        keytype="ENTREZID",
                        columns=c("SYMBOL", "ENTREZID")) %>%
    as_tibble %>% rlang::set_names("entrez_id", "mgi_symbol")
}

m2me <- function (mgisyms) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=mgisyms,
                        keytype="SYMBOL",
                        columns=c("SYMBOL", "ENTREZID")) %>%
    as_tibble %>% rlang::set_names("mgi_symbol", "entrez_id")
}


m2ens <- function (mgisyms) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=mgisyms,
                              keytype="SYMBOL",
                              columns="ENSEMBL")
  df$ENSEMBL
}

ens2m <- function (ensemblids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=ensemblids,
                              keytype="ENSEMBL",
                              columns="SYMBOL")
  df$SYMBOL
}

ens2mens <- function (ensemblids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=ensemblids,
                        keytype="ENSEMBL",
                        columns=c("SYMBOL", "ENSEMBL")) %>%
    as_tibble %>% rlang::set_names("ensembl_id", "mgi_symbol")
}

m2mens <- function (mgisyms) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=mgisyms,
                        keytype="SYMBOL",
                        columns=c("SYMBOL", "ENSEMBL")) %>%
    as_tibble %>% rlang::set_names("mgi_symbol", "ensembl_id")
}


e2ens <- function (entrezids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=entrezids,
                              keytype="ENTREZID",
                              columns="ENSEMBL")
  df$ENSEMBL
}

ens2e <- function (ensemblids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=ensemblids,
                              keytype="ENSEMBL",
                              columns="ENTREZID")
  df$ENTREZID
}

ens2eens <- function (ensemblids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=ensemblids,
                        keytype="ENSEMBL",
                        columns=c("ENTREZID", "ENSEMBL")) %>%
    as_tibble %>% rlang::set_names("ensembl_id", "entrez_id")
}

e2eens <- function (entrezids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=entrezids,
                        keytype="ENTREZID",
                        columns=c("ENTREZID", "ENSEMBL")) %>%
    as_tibble %>% rlang::set_names("entrez_id", "ensembl_id")
}


m2r <- function (mgisyms) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=mgisyms,
                              keytype="SYMBOL",
                              columns="REFSEQ")
  df$REFSEQ
}

r2m <- function (refseqids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=refseqids,
                              keytype="REFSEQ",
                              columns="SYMBOL")
  df$SYMBOL
}

r2mr <- function (refseqids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=refseqids,
                        keytype="REFSEQ",
                        columns=c("SYMBOL", "REFSEQ")) %>%
    as_tibble %>% rlang::set_names("refseq_id", "mgi_symbol")
}

m2mr <- function (mgisyms) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=mgisyms,
                        keytype="SYMBOL",
                        columns=c("SYMBOL", "REFSEQ")) %>%
    as_tibble %>% rlang::set_names("mgi_symbol", "refseq_id")
}



e2r <- function (entrezids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=entrezids,
                              keytype="ENTREZID",
                              columns="REFSEQ")
  df$REFSEQ
}

r2e <- function (refseqids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=refseqids,
                              keytype="REFSEQ",
                              columns="ENTREZID")
  df$ENTREZID
}

r2er <- function (refseqids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=refseqids,
                        keytype="REFSEQ",
                        columns=c("ENTREZID", "REFSEQ")) %>%
    as_tibble %>% rlang::set_names("refseq_id", "entrez_id")
}

e2er <- function (entrezids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=entrezids,
                        keytype="ENTREZID",
                        columns=c("ENTREZID", "REFSEQ")) %>%
    as_tibble %>% rlang::set_names("entrez_id", "refseq_id")
}



ens2r <- function (enstids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=enstids,
                              keytype="ENSEMBLTRANS",
                              columns="REFSEQ")
  df$REFSEQ
}

r2ens <- function (refseqids) {
  df <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=refseqids,
                              keytype="REFSEQ",
                              columns="ENSEMBLTRANS")
  df$ENSEMBLTRANS
}

r2ensr <- function (refseqids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=refseqids,
                        keytype="REFSEQ",
                        columns=c("ENSEMBLTRANS", "REFSEQ")) %>%
    as_tibble %>% rlang::set_names("refseq_id", "enst_id")
}

ens2ensr <- function (enstids) {
  AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=enstids,
                        keytype="ENSEMBLTRANS",
                        columns=c("ENSEMBLTRANS", "REFSEQ")) %>%
    as_tibble %>% rlang::set_names("enst_id", "refseq_id")
}


keepsynhelper <- function (alias, sym) {
  ## if there is exact synonym, return the index for only these (can be many)
  if (any(alias == sym, na.rm=TRUE))
    which(alias == sym)
  ## otherwise, keep everything
  else
    seq_along(alias)
}


## keeps true synonyms in a symbol alias data frame
keepsyn <- function (ttab) {
  do.call(rbind, by(ttab, ttab$ALIAS, function (df) {
    df[keepsynhelper(df$ALIAS, df$SYMBOL), ]
  }))
}



# Stats -------------------------------------------------------------------



my.fisher.test.num <- function(n1a, n1b, n2a, n2b, alternative="t",
                               verbose=FALSE, pvalando=TRUE) {
  fm <- matrix(c(n1a, n1b, n2a, n2b), nrow=2)
  ft <- fisher.test(fm, alternative=alternative)
  if (verbose) {
    return(list(ft=ft, fm=fm, fr=fm[1, ]/(fm[1, ] + fm[2, ])))
  }
  if (pvalando) {
    return(c(`p-value`=ft$p.value, ft$estimate))
  }
  else {
    print(ft$p.value)
  }
}



#' Performs Fisher test for odds difference of condition 1 contingent on
#' condition 2.
#'
#' NAs must be taken into account when formulating the conditions below,
#' if desired.
#'
#' The matrix (returned if \code{verbose=TRUE}) is constructed with cond1
#' and !cond1 for the first and second row, respectively; likewise with cond2
#' and !cond2 for the first and second column, respectively

#'
#' @param df tibble
#' @param cond1 quosure evaluating to a logical for condition 1
#' @param cond2 quosure evaluating to a logical for condition 2
#' @param alternative default="two.sided", see fisher.test().
#' @param ... further arguments to \code{my.fisher.test.num()} such as
#' \code{verbose}
#'
#' @return p value and odds ratio (c1c2/!c1c2) / (c1!c2/!c1!c2)
#' @export
#'
#' @examples
#' df <- data_frame(a=1:10, b=4:13)
#' fisher_test_df(df, quo(a < 4), quo(b < 7))
#' fisher_test_df(df, quo(a > 2), quo(b < 8))
fisher_test_df <- function (df, cond1, cond2, alternative="two.sided",
                            ...) {
  # Contingency table for "outof" version
  #
  #             cond2            | !cond2
  #         -----------------------------
  #   cond1  |  i(df_1, df_2)  | i(df_1, !df_2)
  #   !cond1 |  i(!df_1, df_2) | i(!df_1, !df_2)
  #

  my.fisher.test.num(nrow(dplyr::filter(df, (!! cond1)    & (!! cond2))     ),
                     nrow(dplyr::filter(df, (!(!! cond1)) & (!! cond2))     ),
                     nrow(dplyr::filter(df, (!! cond1)    & (!(!! cond2)) ) ),
                     nrow(dplyr::filter(df, (!(!! cond1)) & (!(!! cond2)) ) ),
                     alternative=alternative, ...
  )

}



#' Perform p value meta analysis using MetaDE
#'
#' @param ... 2 or more vectors of equal lengths of p values
#' @param metameth meta analysis method, defaluts to minP
#' (OR=union-intersection)
#'
#' @return vector of p values
#' @export
#'
#' @examples qmetade(runif(10), runif(10))
qmetade <- function (..., metameth="minP") {
  # metaan <- MetaDE::MetaDE.pvalue(list(p=do.call(cbind, list(...)),
  #                                      bp=NULL),
  #                                 meta.method=metameth)
  # as.numeric(metaan$meta.analysis$pval)
  pvals <- do.call(cbind, list(...))
  if (ncol(pvals) == 1L)
    return(as.vector(pvals))
  else {
    switch(metameth,
           minP = apply(pvals, 1, function (x)
             ifelse(any(is.na(x)), NA, metap::minimump(x)$p)),
           maxP = apply(pvals, 1, function (x)
             ifelse(any(is.na(x)), NA, metap::maximump(x)$p)),
           Fisher = apply(pvals, 1, function (x)
             ifelse(any(is.na(x)), NA, metap::sumlog(x)$p)),
           Stouffer = apply(pvals, 1, function (x)
             ifelse(any(is.na(x)), NA, metap::sumz(x)$p)),
           stop("Incorrect metameth argument"))
  }
}


# Adjusted p values with pre-filtering

#' Compute adjusted p values with pre-filtering
#'
#' @param pvals vector of input p values
#' @param ind logical vector corresponding to pre-filtering condition. NAs are
#'   allowed and reinterprated as FALSE.
#' @param method method selection for p.adjust
#'
#' @return vector of adjusted p values with NAs for non-satisfied pre-filtering
#' @export
#'
#' @examples
#' pvals <- runif(10)
#' effectsize <- runif(10)
#' sel_padj(pvals, effectsize > 0.5)
sel_padj <- function (pvals, ind, method="BH") {
  indnona <- ind & !is.na(ind)
  toret <- rep(NA, length(pvals))
  toret[indnona] <- stats::p.adjust(pvals[indnona], method)
  toret
}



# q values with pre-filtering

#' Compute q values with pre-filtering
#'
#' @param pvals vector of input p values
#' @param ind logical vector corresponding to pre-filtering condition. NAs are
#'   allowed and reinterprated as FALSE. If NULL, all are interpreted as TRUE.
#' @param method method selection for qvalue::pi0est
#' @param ... passed on to qvalue::qvalue (and then to pi0est)
#'
#' @return vector of q values with NAs for non-satisfied pre-filtering
#' @export
#'
#' @examples
#' pvals <- runif(10)
#' effectsize <- runif(10)
#' sel_qval(pvals, effectsize > 0.5)
sel_qval <- function (pvals, ind = NULL, method = "bootstrap", ...) {
  if (is.null(ind))
    ind <- rep(TRUE, length(pvals))
  indnona <- ind & !is.na(ind) & !is.na(pvals)
  toret <- rep(NA, length(pvals))
  if (sum(indnona) > 0L)
    toret[indnona] <- qvalue::qvalue(pvals[indnona],
                                     pi0.method = method,
                                     ...)$qvalues
  toret
}




# proportion.true.disc <- function(qvals) {
# 	qvals <- qvals[order(qvals)]
# 	rho <- (1:length(qvals))/length(qvals)
# 	proportions <- (1 - qvals)*rho
# 	return(list(proportion=max(proportions), qvals.sorted=qvals,
# 							proportions.vec=proportions))
# }



# Circadian reference data ------------------------------------------------

core_clock_genes_noror <- c("Arntl", "Clock", "Npas2",
                            "Per1", "Per2", "Per3",
                            "Cry1", "Cry2",
                            "Nr1d1", "Nr1d2",
                            "Bhlhe40", "Bhlhe41",
                            "Dbp", "Nfil3", "Hlf", "Tef",
                            "Ciart")

core_clock_genes <- c("Arntl", "Clock", "Npas2",
                      "Per1", "Per2", "Per3",
                      "Cry1", "Cry2",
                      "Nr1d1", "Nr1d2", "Rora", "Rorc",
                      "Bhlhe40", "Bhlhe41",
                      "Dbp", "Nfil3", "Hlf", "Tef",
                      "Ciart")



