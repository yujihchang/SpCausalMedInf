#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
})

usage <- function() {
  cat("
Usage:
  Rscript main/main.R \\
    --approach semi|nonpara \\
    --data data.csv \\
    --times 1,3,5,10 \\
    --effect DE|IE|BOTH \\
    --expo 0,1,0,1 \\
    --center center|median|SET \\
    [--setting 0.2,0.5,...] \\
    [--xcols HBSAG,AGE,...]      (semiparametric; optional,其餘欄位將自動當共變數) \\
    [--zacols ZA]                (nonparametric; 指出暴露欄位名，單一欄即可) \\
    [--outdir results]

Input columns:
  必要欄位：A, X1, X2, D。若資料沒有 'time' 欄，請用 --times 指定查表時間。
Notes:
  * approach=semi 會 source('R/YTH_COX_complete v2.r') 並呼叫 YTH_COX_diff(...)
  * approach=nonpara 會 source('R/YTH_haz general as cox.r') 並呼叫 YTH_haz_complete(...)
  * center:
      - center  : 以各共變數平均數代入（傳給函式的 center='TRUE'）
      - median  : 以各共變數中位數代入（在 main 端轉為 center='SET' 並自動給 setting）
      - SET     : 需用 --setting 指定長度 p-1 的數值（p 為 X 的欄數含截距）
\n")
  quit(status = 1)
}

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, default=NULL) {
  hit <- which(args == flag)
  if (length(hit)==0) return(default)
  if (hit == length(args)) return(TRUE)
  val <- args[hit + 1]
  if (startsWith(val, "--")) return(TRUE)
  val
}

if (length(args) == 0 || "--help" %in% args) usage()

approach <- tolower(getArg("--approach"))
if (!(approach %in% c("semi","nonpara"))) usage()

data_path <- getArg("--data"); if (is.null(data_path)) usage()
times_arg <- getArg("--times")
effect_arg <- toupper(getArg("--effect", "BOTH"))
expo_arg   <- getArg("--expo", "0,1,0,1")
center_arg <- tolower(getArg("--center", "center"))
setting_arg <- getArg("--setting", NULL)
xcols_arg  <- getArg("--xcols", NULL)
zacols_arg <- getArg("--zacols", NULL)
outdir <- getArg("--outdir", "results")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 讀資料
dat <- read.csv(data_path, check.names = FALSE)

# time: 若資料內沒有就吃 --times
if (!("time" %in% names(dat))) {
  if (is.null(times_arg)) stop("資料檔沒有 'time' 欄，請用 --times 指定 (逗號分隔)。")
  time <- as.numeric(strsplit(times_arg, ",", fixed=TRUE)[[1]])
} else {
  time <- dat$time
}

# 必要欄位
required <- c("A","X1","X2","D")
miss <- setdiff(required, names(dat))
if (length(miss) > 0) stop(paste0("缺少欄位: ", paste(miss, collapse=", ")))

A  <- dat$A
X1 <- dat$X1
X2 <- dat$X2
D  <- dat$D

# 檢查 X1<=X2
if (any(X1 > X2, na.rm = TRUE)) {
  stop("條件不符：存在 X1 > X2 的觀測值（需滿足 X1 <= X2）。")
}

# 解析 expo/expor
expo_vals <- as.numeric(strsplit(expo_arg, ",", fixed=TRUE)[[1]])
if (length(expo_vals) != 4) stop("--expo 必須是四個數值，例：0,1,0,1 或 1,2,1,2")

# EFFECT
EFFECT <- switch(effect_arg,
  "DE" = "DE",
  "IE" = "IE",
  "BOTH" = c("DE","IE"),
  stop("--effect 僅支援 DE, IE, BOTH")
)

# 來源路徑（依你的檔名）
if (approach == "semi") {
  src_path <- file.path("R", "YTH_COX_complete v2.r")
  if (!file.exists(src_path)) stop("找不到檔案：R/YTH_COX_complete v2.r")
  source(src_path, local = TRUE)

  # 共變數矩陣 X：若指定 --xcols 則用之；否則自動抓除去 (time,A,X1,X2,D) 的其餘欄
  if (!is.null(xcols_arg)) {
    xcols <- strsplit(xcols_arg, ",", fixed=TRUE)[[1]]
    X <- as.matrix(dat[, xcols, drop=FALSE])
  } else {
    covar_cols <- setdiff(names(dat), c("time","A","X1","X2","D"))
    if (length(covar_cols) == 0) {
      X <- matrix(1, nrow(dat), 1)
      colnames(X) <- "(Intercept)"
    } else {
      X <- as.matrix(cbind(1, dat[, covar_cols, drop=FALSE]))
      colnames(X)[1] <- "(Intercept)"
    }
  }
  p <- ncol(X)

  # center 轉換成函式可接受的參數
  # 你的函式支援 center='TRUE' (用平均)、'FALSE'、'SET'
  center_flag <- switch(center_arg,
    "center" = "TRUE",
    "median" = "SET",
    "set"    = "SET",
    stop("--center 僅支援 center|median|SET")
  )

  if (center_flag == "SET") {
    if (tolower(center_arg) == "median") {
      # 自動以中位數做 setting（僅非截距欄）
      if (p > 1) {
        setting <- apply(as.matrix(X[,-1,drop=FALSE]), 2, median, na.rm = TRUE)
      } else {
        setting <- numeric(0)
      }
    } else {
      if (is.null(setting_arg)) stop("center='SET' 時，請用 --setting 提供 p-1 個值（逗號分隔）")
      setting <- as.numeric(strsplit(setting_arg, ",", fixed=TRUE)[[1]])
    }
  } else {
    setting <- NULL
  }

  # 呼叫半參數法
  res <- YTH_COX_diff(
    time = time,
    A    = A,
    X1   = X1,
    X2   = X2,
    D    = D,
    X    = X,
    expor  = expo_vals,
    center = center_flag,
    setting = setting,
    EFFECT  = EFFECT,
    saving  = "FALSE"
  )

  out <- data.frame(t = time)
  if ("DE" %in% EFFECT) { out$DE <- res$DE; out$DE_sd <- res$DE.sd }
  if ("IE" %in% EFFECT) { out$IE <- res$IE; out$IE_sd <- res$IE.sd }
  out$LA_Sz0_za2 <- res$LA_Sz0_za2
  out$LA_Sz1_za2 <- res$LA_Sz1_za2
  out$LA_Sz0_za1 <- res$LA_Sz0_za1
  out$LA_Sz1_za1 <- res$LA_Sz1_za1

  outfile <- file.path(outdir, sprintf("semi_YTH_COX_diff_%s.csv", tolower(effect_arg)))
  write.csv(out, outfile, row.names = FALSE)
  cat(sprintf("Semiparametric results => %s\n", outfile))

} else {
  # Nonparametric
  src_path <- file.path("R", "YTH_haz general as cox.r")
  if (!file.exists(src_path)) stop("找不到檔案：R/YTH_haz general as cox.r")
  source(src_path, local = TRUE)

  if (is.null(zacols_arg)) stop("nonparametric 需要 --zacols 指出暴露欄位名 (單一欄)")
  Za <- dat[[zacols_arg]]

  res <- YTH_haz_complete(
    t = time,
    A = A,
    X1 = X1,
    X2 = X2,
    D = D,
    Za = Za,
    expo = expo_vals,
    EFFECT = EFFECT
  )

  out <- data.frame(t = time)
  if ("DE" %in% EFFECT) { out$DE <- res$DE; out$DE_sd <- res$DE.sd }
  if ("IE" %in% EFFECT) { out$IE <- res$IE; out$IE_sd <- res$IE.sd }

  outfile <- file.path(outdir, sprintf("nonpara_YTH_haz_complete_%s.csv", tolower(effect_arg)))
  write.csv(out, outfile, row.names = FALSE)
  cat(sprintf("Nonparametric results => %s\n", outfile))
}
