source('../MDICC/NetworkFusion.R')
library(Rcpp)
library(parallel)
library(Matrix)

setwd("../MDICC/")
# system("R CMD SHLIB projsplx_R.c")
dyn.load("projsplx_R.dll")

# connect anaconda environment
library(reticulate)
use_virtualenv("base")
# the path of python 
Sys.setenv(RETICULATE_PYTHON="C:\\Users\\*\\anaconda3\\python.exe") 
use_python("C:\\Users\\*\\anaconda3\\python.exe")
py_config()
py_available()
source_python("LocalAffinityMatrix.py")
source_python("score.py")
source_python("label.py")

# read data
setwd("../MDICC/data/brca/") # data path
list <- list.files()
data <- data.frame()
data1 <- list()
X <- list()
for(i in list){
  path <- i
  data <- read.csv(file = path, header = TRUE)
  # rownames(data) <- data$X
  data <- data[-1]
  data11 <- as.matrix(data)
  data1[[i]] = scale(data11, center=TRUE, scale=TRUE) 
  data2 = t(data1[[i]])
  d1 = dist(data2)
  d1 = as.matrix(d1)
  X[[i]] <- d1
}


# parameter setting
k1 = 18 # the neighbor of affinity matrix
k2 = 42 # 
k3 = 2  # number of cluster
c  = 3  # c = k3(c>2) or c = 3(c=2)

# calculate affinity matrix
aff = list()
for(i in 1:3){
  a = as.matrix(X[[i]])
  xxx = testaff(a,k1)
  aff[[i]] = xxx
}

# network fusion
test = MDICC(aff,c = c,k = k2)
test_S = as.matrix(test)

# result
# label path 
score = MDICCscore(test_S,k3,'.../MDICC/data/label.csv','label1')
names(score) = c('RI','ARI','NMI','Accu','F1')
label = MDICClabel(test_S,k3)
MDICCresult = list()
MDICCresult[['score']] = score
MDICCresult[['label']] = label


