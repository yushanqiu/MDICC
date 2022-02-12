Package: MDICC
Type: Package
Title: MDICC: Novel method for multi-omics data integration and cancer subtype identification
Version: 1.0
Author: Yang Ying <m13265355293@163.com>
Maintainer: Yushan Qiu<yushan.qiu@szu.edu.cn>
Description: This pakege implemented the MDICC algorithm that fusion different omics data, which aims to identify cancer subtypes.
        
### Environment: 
Anaconda, R>= 4.1.0;
### R package: 
Rcpp;   parallel;  Matrix; reticulate >= 1.22;
### Introduction: 
The algorithm is implemented by Python and R.  You need to install the Anaconda environment and Rstudio, and the R version used in this experiment is 4.1.0.  The source code is mainly written in R.  The experimental data are BRCA omics data and stored in the folder "data".  The steps are as follows：

###### Import MIDCC algorithms from data sources.  
> source('../MDICC/NetworkFusion.R')   
###### Import the projsplx_r.dll file.
> setwd("../MDICC/")    
> dyn.load("projsplx_R.dll")  
##### If that fails, you can run the following code.
> system("R CMD SHLIB projsplx_R.c")  
##### Connect Anaconda environment.
> library(reticulate)  
> use_virtualenv("base")  
> Sys.setenv(RETICULATE_PYTHON="C:\\Users\\*\\anaconda3\\python.exe")   
> use_python("C:\\Users\\*\\anaconda3\\python.exe")  
> py_config()  
> py_available()  
##### If the command is successfully executed, TRUE is displayed. Then import the Python files.
> source_python("LocalAffinityMatrix.py")  
> source_python("score.py")  
> source_python("label.py") 
##### Where, score.py is the file that calculates ARI and NMI, and label.py is the file that gets the label after clustering. Next, you can perform the fusion analysis of multi-omics data according to the demo file.


 - If you have any problem, please contact Ying Yang (m13265355293@163.com) and Yushan Qiu (yushan.qiu@szu.edu.cn)!



