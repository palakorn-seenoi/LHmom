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



#------------ การส่ง R Package ให้ อ.พัค ------------

#วิธีที่ 1: บีบอัดเป็นไฟล์ .zip (รวดเร็วและง่ายที่สุด)
devtools::install_local("ที่อยู่ไฟล์/LHmom.zip")

#วิธีที่ 2: สร้างเป็นไฟล์แพ็กเกจมาตรฐาน (.tar.gz)
devtools::build()
install.packages("ที่อยู่ไฟล์/LHmom_0.1.0.tar.gz", repos = NULL, type = "source")

#วิธีที่ 3: อัปโหลดขึ้น GitHub (มืออาชีพและอัปเดตง่ายที่สุด)
devtools::install_github("ชื่อ_username_ของคุณ/LHmom")
devtools::install_github("palakorn-seenoi/LHmom")


# use_r -------------------------------------------------------------------

usethis::use_r("lh.pargno.R")

devtools::document()

# use_test ----------------------------------------------------------------

#usethis::use_test("lh.pargno.test.R")

#devtools::document()

# use_import_from ----------------------------------------------------------------

#use_package("lmomco")
#use_package("nleqslv")

use_import_from("lmomco", "pargno")
use_import_from("nleqslv", "nleqslv")

devtools::document()

#load_all()
#rm(list = c("lh.pargno", "lhmom.gno", "lhmoms"))
#-----------------------------------------------

#use_description()

use_package_doc()

devtools::test()
devtools::check()
load_all()
devtools::build()
devtools::build(binary=TRUE)

#----------- Checking Package ---------------------------------------

library(LHmom)
hello()
data(package = "LHmom")
ls("package:LHmom")
data("rain_khonkaen")

help(package = "LHmom")

Imax
lh.pargno
lhmom.gno
lhmoms
pargno


#----------- sample Running ---------------------------------------

data("rain_khonkaen")
rain_khonkaen
data <- rain_khonkaen$max_rain
data

eta=0
pargno(lmoms(data))
lh.pargno(data,eta=eta, opt=FALSE)
lh.pargno(data,eta=eta, opt=TRUE)

for (eta in 0:4){
  lhold= lh.pargno(data,eta=eta, opt=FALSE)
  lhnew= lh.pargno(data,eta=eta, opt=TRUE)
  cat("eta, para.old =", eta, lhold$para,"\n")
  cat("eta, para.new =", eta, lhnew$para,"\n")
}


