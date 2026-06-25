library(LHmom)

data(bangkok1)

# Calculate sample LH-moments at eta = 1
sample_lh <- lhmoms(bangkok1$rainfall, eta = 1)
print(sample_lh$lambdas)
##        lhmom-1  lhmom-2  lhmom-3  lhmom-4  lhmom-5
## eta=1 126.5447 18.86068 5.220361 3.037946 1.666922

# Estimate GEV parameters (eta = 1)
fit_gev <- lh.pargev(bangkok1$rainfall, eta = 1)
print(fit_gev$para)
#>          xi       alpha           k
#> 87.64384852 28.23812277 -0.06139394

# Calculate theoretical LH-moments from the fitted parameters
theo_lh_gev <- lhmom.gev(fit_gev$para, eta = 1)
print(theo_lh_gev$lambdas)
#>    LHmom-1    LHmom-2    LHmom-3    LHmom-4
#> 126.544685  18.860680   5.229035   3.131879



#--------------------------------------------------------------------------------------

# Generate Q-Q plot for GEV ( eta = 1) with 95% bootstrap CI
lh.qqplot( bangkok1$rainfall, fit_gev, main = "GEV Q-Q Plot (eta = 1)",
           ci = TRUE, ci_level = 0.95)
#--------------------------------------------------------------------------------------


# Perform Wang's goodness-of-fit test for the GEV distribution
gof_test <- wang.test.lhgev(bangkok1$rainfall)
print(gof_test)
##    eta       z.test cond.sigma   p.value
##  1   0 -0.149710710 0.03581892 0.8809929
##  2   1 -0.167327234 0.02976428 0.8671126
##  3   2 -0.170480325 0.02714540 0.8646324
##  4   3  0.006399968 0.02572831 0.9948936
##  5   4  0.170936924 0.02478531 0.8642734


#--------------------------------------------------------------------------------------



# 1. เตรียมข้อมูล
data(bangkok1)
x <- bangkok1$rainfall
x

# 2. ตั้งค่าหน้าจอสำหรับ 5 กราฟ (2 แถว 3 คอลัมน์)
pdf("Figure_GEV_All_Etas.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))

# 3. กำหนดค่า eta และวนลูป
etas <- 0:4
fits_gev <- list()

for (e in etas) {
  # ฟิตพารามิเตอร์ GEV
  fits_gev[[as.character(e)]] <- lh.pargev(x, eta = e)

  # กำหนดชื่อกราฟ
  plot_title <- paste("GEV Q-Q Plot (eta =", e, ")")

  # วาด Q-Q Plot พร้อมแถบ CI 95% และแสดงค่า RMSE
  lh.qqplot(x, fits_gev[[as.character(e)]],
            main = plot_title,
            ci = TRUE, ci_level = 0.95)
}
dev.off()





#---------------------------------------------------------------------
