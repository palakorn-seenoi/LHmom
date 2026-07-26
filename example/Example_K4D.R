library(LHmom)

data(bangkok1)

# Calculate sample LH-moments at eta = 2
sample_lh <- lhmoms(bangkok1$rainfall, eta = 2)
print(sample_lh$lambdas)
##         lhmom-1  lhmom-2  lhmom-3  lhmom-4  lhmom-5
##  eta=2 139.1185 18.22076 5.464799 2.951366 1.867736
print(sample_lh$ratios)
##          lht-1    lht-2     lht-3     lht-4     lht-5
##  eta=2      NA 0.130973 0.2999215 0.1619782 0.1025059

# Estimate Kappa parameters using sample LH-moments (eta = 2)
fit_kap <- lh.parkap(bangkok1$rainfall, eta = 2)
print(fit_kap$para)
##          mu      sigma          k          h
##  84.4133018 30.7063710 -0.0312735  0.1986932



# Calculate theoretical LH-moments from the fitted parameters
theo_lh_kap <- lhmom.kap(fit_kap$para, eta = 2)
print(theo_lh_kap$lambdas)
##    LHmom-1    LHmom-2    LHmom-3    LHmom-4
## 139.118472  18.220761   5.464818   2.951380
print(theo_lh_kap$ratios)
##   LHtau-1   LHtau-2   LHtau-3   LHtau-4
##        NA 0.1309730 0.2999226 0.1619790


print(fit_kap$ifailtext)

print(fit_kap$precision)


#--------------------------------------------------------------------------------------

# Generate Q-Q plot for Kappa (eta = 2) with 95% bootstrap CI
lh.qqplot( bangkok1$rainfall, fit_kap,
           main = "Kappa Q-Q Plot (eta = 2)",
           ci = TRUE, ci.level = 0.95)
#--------------------------------------------------------------------------------------


# 1. เตรียมข้อมูล
data(bangkok1)
x <- bangkok1$rainfall
x

# 2. ตั้งค่าหน้าจอสำหรับ 5 กราฟ (2 แถว 3 คอลัมน์ เหมือน GEV)
pdf("Figure_Kappa_All_Etas.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))

# 3. กำหนดค่า eta และวนลูป
etas_kap <- 0:4
fits_kap <- list()

for (e in etas_kap) {
  # ฟิตพารามิเตอร์ 4-parameter Kappa
  fits_kap[[as.character(e)]] <- lh.parkap(x, eta = e)

  # กำหนดชื่อกราฟ
  plot_title <- paste("Kappa Q-Q Plot (eta =", e, ")")

  # วาด Q-Q Plot พร้อมแถบ CI 95% (เอา RMSE ออกแล้ว)
  lh.qqplot(x, fits_kap[[as.character(e)]],
            main = plot_title,
            ci = TRUE, ci.level = 0.95)
}
dev.off()





