install.packages(c("devtools", "roxygen2", "usethis"))
install.packages("devtools")
library(devtools)
library(roxygen2)
library(usethis)


usethis::create_package("D:/SoftwareX/LHmom/GitHub/LHmom")


# คำสั่งนี้จะสร้างไฟล์คู่มือ (Documentation) อัตโนมัติจากคอมเมนต์ในข้อ 3
devtools::document()

# คำสั่งนี้จะทำการติดตั้ง Library นี้ลงในเครื่องของคุณ
devtools::install()

#devtools::load_all() หรือใช้คีย์ลัด Ctrl + Shift + L
#ระบบจะจำลองการโหลดโค้ดทั้งหมดในโฟลเดอร์ R/ ขึ้นมาให้ใช้งานได้ทันทีโดยไม่ต้องติดตั้งใหม่
devtools::load_all()





# use_r -------------------------------------------------------------------

usethis::use_r("lhmoms")
usethis::use_r("initk.lhme.R")  # Do not to show this function
usethis::use_r("lh.pargev.R")
usethis::use_r("lh.pargno.R")
usethis::use_r("lh.pargpa_glo.R")
usethis::use_r("lh.pargum.R")
usethis::use_r("lh.park3d.hfix.R")
usethis::use_r("lh.park3d.kfix.R")
usethis::use_r("lh.park4d.R")
usethis::use_r("lh.parpe3.R")
usethis::use_r("wang.test.lhgev.R")
usethis::use_r("theo.lhmom.others.R")
usethis::use_r("datasets.R")

devtools::document()

# use_test ----------------------------------------------------------------

#usethis::use_test("lh.pargno.test.R")

#devtools::document()

# use_import_from ----------------------------------------------------------------

usethis::use_testthat()

usethis::use_package("lmomco")
usethis::use_package("nleqslv")
usethis::use_package("stats")


use_import_from("lmomco", "lmoms")
use_import_from("lmomco", "lmomkap")

use_import_from("lmomco", "pargev")
use_import_from("lmomco", "pargno")
use_import_from("lmomco", "parkap")
use_import_from("lmomco", "parpe3")
use_import_from("lmomco", "vec2par")

use_import_from("nleqslv", "nleqslv")

devtools::document()

#load_all()
#rm(list = c("lh.pargno", "lhmom.gno", "lhmoms"))
#-----------------------------------------------

usethis::use_build_ignore("Setting_LHmom.R")

use_description()

use_package_doc()

devtools::test()
devtools::check()
load_all()
devtools::build()
devtools::build(binary=TRUE)



#----------- Checking Package ---------------------------------------

library(LHmom)
data(package = "LHmom")
ls("package:LHmom")
data("rain_khonkaen")

help(package = "LHmom")
?lhmom.k3d.hfix
?lhmom.gev
?lh.pargev
?lh.pargno
?lhmom.gno
?Imax      #internal
?lh.parglo
?lh.pargpa
?lh.pargum
?lh.park3d.hfix
?lh.parggd
?lh.park3d.kfix
?lh.parkap
?lhmom.kap
?initkh.k4d #internal
?lh.parpe3
?lhmom.pe3
?cal_Omega #internal
?cal_Lam   #internal
?wang.test.lhgev
?lhmom.glo
?lhmom.gpa
?lhmom.ggd
?lhmom.k3d.hfix
?lhmom.k3d.kfix
?lhmom.gum

lmoms(data)
lhmoms(data, eta = 1)
pargno


#----------- Sample Running ---------------------------------------

#--------------- data 1 ---------------
data("rain_khonkaen")
rain_khonkaen
data <- rain_khonkaen$max_rain
data

eta=0
pargno(lmoms(data))
lh.pargno(data, eta=eta, opt=FALSE)
lh.pargno(data, eta=eta, opt=TRUE)



for (eta in 0:4){
  lhold= lh.pargno(data, eta=eta, opt=FALSE)
  lhnew= lh.pargno(data, eta=eta, opt=TRUE)
  cat("eta, para.old =", eta, lhold$para,"\n")
  cat("eta, para.new =", eta, lhnew$para,"\n")
}

#--------------- data 2 ---------------
data("rain_bangkok")
rain_bangkok
data <- rain_bangkok$max_rain
data


#---------- Sample Running ---------------------------------------

eta=1
####### Theoretical LH-moments #######
#1. GEV : Generalized extreme value distribution (gev)
lhmom.gev(, eta = eta)
#2. GGD : generalized Gumbel distribution (ggd)
#3. GLO : Generalized logistic distribution (glo)
#4. GNO : generalized normal distribution (gno)
lhmom.gno(, eta = eta)
#5. GPA : Generalized Pareto distribution (gpa)
#6. GUM : Gumbel distribution (gum)
#7. k3d.hfix : Three-parameter kappa distribution with h-fixed (h3d_hfix)
#8. k3d.kfix : Three-parameter kappa distribution with k-fixed (k3d_kfix)
#9. kap : Four-parameter kappa distribution (kap)
lhmom.kap(, eta = eta)
#10. pe3 : Pearson type-3 distribution (pe3)
lhmom.pe3(, eta = eta)


####### Parameter Estimation LH-moments #######
eta=1
#1. GEV : Generalized extreme value distribution (gev)
lh.pargev(data, eta = eta)
#2. GGD : generalized Gumbel distribution (ggd)
lh.parggd(data, eta = eta)
#3. GLO : Generalized logistic distribution (glo)
lh.parglo(data, eta = eta)
#4. GNO : generalized normal distribution (gno)
lh.pargno(data, eta = eta)
#5. GPA : Generalized Pareto distribution (gpa)
lh.pargpa(data, eta = eta)
#6. GUM : Gumbel distribution (gum)
lh.pargum(data, eta = eta)
#7. k3d.hfix : Three-parameter kappa distribution with h-fixed (h3d_hfix)
lh.park3d.hfix(data, eta = 0)
lh.park3d.hfix(data, eta = 1)
#8. k3d.kfix : Three-parameter kappa distribution with k-fixed (k3d_kfix)
lh.park3d.kfix(data, eta = eta)
#9. kap : Four-parameter kappa distribution (kap)
lh.parkap(data, eta = 0)
lh.parkap(data, eta = 1)
#10. pe3 : Pearson type-3 distribution (pe3)
lh.parpe3(data, eta = 2)



##> lh.parkap(data, eta = eta)
##Error in park3h(lmoms(data), hfix = -2) :
##  could not find function "park3h"



####### Sample LH-moments #######
lhmoms(data, eta = eta)


####### wang.test.lhgev #######
wang.test.lhgev(data)


initk
initkh.k4d


#------------ การส่ง R Package ให้ อ.พัค ------------

#วิธีที่ 1: บีบอัดเป็นไฟล์ .zip (รวดเร็วและง่ายที่สุด)
devtools::install_local("ที่อยู่ไฟล์/LHmom.zip")

#วิธีที่ 2: สร้างเป็นไฟล์แพ็กเกจมาตรฐาน (.tar.gz)
devtools::build()
install.packages("ที่อยู่ไฟล์/LHmom_0.1.0.tar.gz", repos = NULL, type = "source")

#วิธีที่ 3: อัปโหลดขึ้น GitHub (มืออาชีพและอัปเดตง่ายที่สุด)
devtools::install_github("ชื่อ_username_ของคุณ/LHmom")
devtools::install_github("palakorn-seenoi/LHmom")




#------------ Manual of Package ------------

devtools::document()   # man/

#install.packages("tinytex") # install LaTex
tinytex::install_tinytex()
Sys.which("pdflatex")
devtools::build_manual() # PDF

#install.packages("pkgdown") # HTML
usethis::use_pkgdown()

pkgdown::build_site() # ระบบจะสร้างโฟลเดอร์ใหม่ขึ้นมาในโปรเจกต์ชื่อว่า docs/ ซึ่งข้างในจะเต็มไปด้วยไฟล์หน้าเว็บ .html
