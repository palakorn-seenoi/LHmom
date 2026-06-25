library(LHmom)
library(lmomco)
data(bangkok)

# Calculate sample LH-moments at eta = 1
sample_lh <- lhmoms(bangkok$rainfall, eta = 1, nmom = 3)
print(sample_lh$lambdas)

# Estimate GEV parameters using LH-moments (eta = 1)
fit_gev <- lh.pargev(bangkok$rainfall, eta = 1)
print(fit_gev$para)

theo_lh <- lhmom.gev(fit_gev$para, eta = 1)
print(theo_lh$lambdas)

# Perform Wang's goodness-of-fit test for the GEV distribution
gof_test <- wang.test.lhgev(bangkok$rainfall)
print(gof_test)


lh.parggd(bangkok$rainfall)
lh.pargum(bangkok$rainfall)
wang.test.lhgev(bangkok$rainfall)
#---------------------------------------------------------------
etas <- 0:4

# 1. Sample LH-moments (lhmoms)
sample_df <- do.call(rbind, lapply(etas, function(e) {
  res <- lhmoms(bangkok$rainfall, eta = e)
  data.frame(eta = e,
             L1 = res$lambdas[1], L2 = res$lambdas[2],
             L3 = res$lambdas[3], L4 = res$lambdas[4])
}))
print(sample_df)
##   eta       L1       L2       L3       L4
## 1   0 71.59733 13.56586 3.125191 2.457813
## 2   1 85.16320 12.87224 2.910361 2.012546
## ...

# 2. Estimate GEV parameters (lh.pargev)
para_df <- do.call(rbind, lapply(etas, function(e) {
  res <- lh.pargev(bangkok$rainfall, eta = e)
  data.frame(eta = e,
             xi = res$para[1], alpha = res$para[2], k = res$para[3])
}))
print(para_df)
##   eta       xi    alpha           k
## 1   0 62.14384 19.23812 -0.10139394
## 2   1 64.64384 18.23812 -0.06139394
## ...

# 3. Theoretical LH-moments (lhmom.gev)
theo_df <- do.call(rbind, lapply(etas, function(e) {
  para <- lh.pargev(bangkok$rainfall, eta = e)$para
  res <- lhmom.gev(para, eta = e)
  data.frame(eta = e,
             theo_L1 = res$lambdas[1], theo_L2 = res$lambdas[2],
             theo_L3 = res$lambdas[3], theo_L4 = res$lambdas[4])
}))
print(theo_df)
##   eta  theo_L1  theo_L2  theo_L3  theo_L4
## 1   0 71.59733 13.56586 3.125191 2.315243
## 2   1 85.16320 12.87224 2.910361 1.915243
## ...



#---------------------------------------







# 1. โหลดชุดข้อมูล
data(bangkok)
x <- bangkok$rainfall

# 2. ประมาณค่าพารามิเตอร์ GEV ด้วย eta ต่างๆ
fit_eta0 <- lh.pargev(x, eta = 0) # L-moments ปกติ
fit_eta1 <- lh.pargev(x, eta = 1) # LH-moments ระดับ 1
fit_eta2 <- lh.pargev(x, eta = 2) # LH-moments ระดับ 2
fit_eta3 <- lh.pargev(x, eta = 3) # LH-moments ระดับ 2
fit_eta4 <- lh.pargev(x, eta = 4) # LH-moments ระดับ 2

# 3. เตรียมช่วงข้อมูลแกน X สำหรับวาดกราฟ (ครอบคลุมตั้งแต่ค่าน้อยสุดถึงมากสุด)
x_seq <- seq(min(x) * 0.8, max(x) * 1.2, length.out = 200)

# คำนวณค่า Probability Density Function (PDF) ด้วยฟังก์ชันจาก lmomco
# (เนื่องจาก LHmom มีโครงสร้าง Output ตรงกับ lmomco จึงใช้งานร่วมกันได้ทันที)
pdf_eta0 <- pdfgev(x_seq, fit_eta0)
pdf_eta1 <- pdfgev(x_seq, fit_eta1)
pdf_eta2 <- pdfgev(x_seq, fit_eta2)
pdf_eta3 <- pdfgev(x_seq, fit_eta2)
pdf_eta4 <- pdfgev(x_seq, fit_eta2)

# 4. สร้างและบันทึกกราฟเป็นไฟล์ PDF สไตล์วิชาการ
#pdf("figure_gev_density.pdf", width = 8, height = 6)

# วาดกราฟ Histogram ของข้อมูลจริงเป็นพื้นหลัง
hist(x, probability = TRUE, col = "gray90", border = "white",
     main = "GEV Density Comparison for Bangkok Rainfall",
     xlab = "Annual Maximum Rainfall (mm)", ylab = "Density",
     ylim = c(0, max(pdf_eta0, pdf_eta1, pdf_eta2) * 1.1))

# เพิ่มเส้นโค้ง Density ของแต่ละ eta
lines(x_seq, pdf_eta0, col = "blue", lwd = 2, lty = 1)
lines(x_seq, pdf_eta1, col = "red", lwd = 2, lty = 2)
lines(x_seq, pdf_eta2, col = "darkgreen", lwd = 2, lty = 3)
lines(x_seq, pdf_eta3, col = "brown", lwd = 2, lty = 4)
lines(x_seq, pdf_eta4, col = "black", lwd = 2, lty = 5)

# เพิ่ม Legend อธิบายเส้น
legend("topright",
       legend = c(expression(eta == 0~"(Ordinary L-moments)"),
                  expression(eta == 1~"(LH-moments)"),
                  expression(eta == 2~"(LH-moments)"),
                  expression(eta == 3~"(LH-moments)"),
                  expression(eta == 4~"(LH-moments)")),
       col = c("blue", "red", "darkgreen","brown","black"),
       lwd = 2, lty = c(1, 2, 3, 4, 5), bty = "n")

#dev.off() # ปิดคำสั่งและบันทึกไฟล์






#---------------------------------------------------------------
# khonkaen
library(LHmom)
data("khonkaen")
khonkaen

# Calculate sample LH-moments at eta = 1
sample_lh <- lhmoms(khonkaen$rainfall, eta = 1)
print(sample_lh$lambdas)

# Estimate GEV parameters using LH-moments (eta = 1)
fit_gev <- lh.pargev(khonkaen$rainfall, eta = 1)
print(fit_gev$para)

theo_lh <- lhmom.gev(fit_gev$para, eta = 1)
print(theo_lh$lambdas)

# Perform Wang's goodness-of-fit test for the GEV distribution
gof_test <- wang.test.lhgev(khonkaen$rainfall)
print(gof_test)


#---------------------------------------


# 1. โหลดชุดข้อมูล
data(khonkaen)
x <- khonkaen$rainfall

# 2. ประมาณค่าพารามิเตอร์ GEV ด้วย eta ต่างๆ
fit_eta0 <- lh.pargev(x, eta = 0) # L-moments ปกติ
fit_eta1 <- lh.pargev(x, eta = 1) # LH-moments ระดับ 1
fit_eta2 <- lh.pargev(x, eta = 2) # LH-moments ระดับ 2
fit_eta3 <- lh.pargev(x, eta = 3) # LH-moments ระดับ 2
fit_eta4 <- lh.pargev(x, eta = 4) # LH-moments ระดับ 2

# 3. เตรียมช่วงข้อมูลแกน X สำหรับวาดกราฟ (ครอบคลุมตั้งแต่ค่าน้อยสุดถึงมากสุด)
x_seq <- seq(min(x) * 0.8, max(x) * 1.2, length.out = 200)

# คำนวณค่า Probability Density Function (PDF) ด้วยฟังก์ชันจาก lmomco
# (เนื่องจาก LHmom มีโครงสร้าง Output ตรงกับ lmomco จึงใช้งานร่วมกันได้ทันที)
pdf_eta0 <- pdfgev(x_seq, fit_eta0)
pdf_eta1 <- pdfgev(x_seq, fit_eta1)
pdf_eta2 <- pdfgev(x_seq, fit_eta2)
pdf_eta3 <- pdfgev(x_seq, fit_eta2)
pdf_eta4 <- pdfgev(x_seq, fit_eta2)

# 4. สร้างและบันทึกกราฟเป็นไฟล์ PDF สไตล์วิชาการ
#pdf("figure_gev_density.pdf", width = 8, height = 6)

# วาดกราฟ Histogram ของข้อมูลจริงเป็นพื้นหลัง
hist(x, probability = TRUE, col = "gray90", border = "white",
     main = "GEV Density Comparison for Khon Kaen Rainfall",
     xlab = "Annual Maximum Rainfall (mm)", ylab = "Density",
     ylim = c(0, max(pdf_eta0, pdf_eta1, pdf_eta2) * 1.1))

# เพิ่มเส้นโค้ง Density ของแต่ละ eta
lines(x_seq, pdf_eta0, col = "blue", lwd = 2, lty = 1)
lines(x_seq, pdf_eta1, col = "red", lwd = 2, lty = 2)
lines(x_seq, pdf_eta2, col = "darkgreen", lwd = 2, lty = 3)
lines(x_seq, pdf_eta3, col = "brown", lwd = 2, lty = 4)
lines(x_seq, pdf_eta4, col = "black", lwd = 2, lty = 5)

# เพิ่ม Legend อธิบายเส้น
legend("topright",
       legend = c(expression(eta == 0~"(Ordinary L-moments)"),
                  expression(eta == 1~"(LH-moments)"),
                  expression(eta == 2~"(LH-moments)"),
                  expression(eta == 3~"(LH-moments)"),
                  expression(eta == 4~"(LH-moments)")),
       col = c("blue", "red", "darkgreen","brown","black"),
       lwd = 2, lty = c(1, 2, 3, 4, 5), bty = "n")

#dev.off() # ปิดคำสั่งและบันทึกไฟล์



#---------------------------------------
#sarakham

# 1. โหลดชุดข้อมูล
data("sarakham")
sarakham$temperature

# Calculate sample LH-moments at eta = 1
sample_lh <- lhmoms(sarakham$temperature, eta = 1)
print(sample_lh$lambdas)

# Estimate GEV parameters using LH-moments (eta = 1)
fit_gev <- lh.pargev(sarakham$temperature, eta = 1)
print(fit_gev$para)

theo_lh <- lhmom.gev(fit_gev$para, eta = 1)
print(theo_lh$lambdas)

# Perform Wang's goodness-of-fit test for the GEV distribution
gof_test <- wang.test.lhgev(sarakham$temperature)
print(gof_test)


#---------------------------------------


# 1. โหลดชุดข้อมูล
data("sarakham")
x <- sarakham$temperature

# 2. ประมาณค่าพารามิเตอร์ GEV ด้วย eta ต่างๆ
fit_eta0 <- lh.pargev(x, eta = 0) # L-moments ปกติ
fit_eta1 <- lh.pargev(x, eta = 1) # LH-moments ระดับ 1
fit_eta2 <- lh.pargev(x, eta = 2) # LH-moments ระดับ 2
fit_eta3 <- lh.pargev(x, eta = 3) # LH-moments ระดับ 2
fit_eta4 <- lh.pargev(x, eta = 4) # LH-moments ระดับ 2

# 3. เตรียมช่วงข้อมูลแกน X สำหรับวาดกราฟ (ครอบคลุมตั้งแต่ค่าน้อยสุดถึงมากสุด)
x_seq <- seq(min(x) * 0.8, max(x) * 1.2, length.out = 200)

# คำนวณค่า Probability Density Function (PDF) ด้วยฟังก์ชันจาก lmomco
# (เนื่องจาก LHmom มีโครงสร้าง Output ตรงกับ lmomco จึงใช้งานร่วมกันได้ทันที)
pdf_eta0 <- pdfgev(x_seq, fit_eta0)
pdf_eta1 <- pdfgev(x_seq, fit_eta1)
pdf_eta2 <- pdfgev(x_seq, fit_eta2)
pdf_eta3 <- pdfgev(x_seq, fit_eta2)
pdf_eta4 <- pdfgev(x_seq, fit_eta2)

# 4. สร้างและบันทึกกราฟเป็นไฟล์ PDF สไตล์วิชาการ
#pdf("figure_gev_density.pdf", width = 8, height = 6)

# วาดกราฟ Histogram ของข้อมูลจริงเป็นพื้นหลัง
hist(x, probability = TRUE, col = "gray90", border = "white",
     main = "GEV Density Comparison for Maha Sarakham Temerature",
     xlab = "Annual Maximum Rainfall (mm)", ylab = "Density",
     ylim = c(0, max(pdf_eta0, pdf_eta1, pdf_eta2) * 1.1))

# เพิ่มเส้นโค้ง Density ของแต่ละ eta
lines(x_seq, pdf_eta0, col = "blue", lwd = 2, lty = 1)
lines(x_seq, pdf_eta1, col = "red", lwd = 2, lty = 2)
lines(x_seq, pdf_eta2, col = "darkgreen", lwd = 2, lty = 3)
lines(x_seq, pdf_eta3, col = "brown", lwd = 2, lty = 4)
lines(x_seq, pdf_eta4, col = "black", lwd = 2, lty = 5)

# เพิ่ม Legend อธิบายเส้น
legend("topright",
       legend = c(expression(eta == 0~"(Ordinary L-moments)"),
                  expression(eta == 1~"(LH-moments)"),
                  expression(eta == 2~"(LH-moments)"),
                  expression(eta == 3~"(LH-moments)"),
                  expression(eta == 4~"(LH-moments)")),
       col = c("blue", "red", "darkgreen","brown","black"),
       lwd = 2, lty = c(1, 2, 3, 4, 5), bty = "n")

#dev.off() # ปิดคำสั่งและบันทึกไฟล์


