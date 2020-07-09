#' check the environment for necessary packages, and install the missing packages
#'
#' @examples
#' env_check()

env_check <- function(){
  if("mgcv" %in% rownames(installed.packages()) == FALSE)
  {install.packages("mgcv")
  }
  library("mgcv")

  if("MASS" %in% rownames(installed.packages()) == FALSE)
  {install.packages("MASS")
  }
  library("MASS")

  if("gtools" %in% rownames(installed.packages()) == FALSE)

  {install.packages("gtools")
  }
  library("gtools")

  if("glmnet" %in% rownames(installed.packages()) == FALSE)
  {install.packages("glmnet")
  }
  library("glmnet")

  if("Matrix" %in% rownames(installed.packages()) == FALSE)
  {install.packages("Matrix")
  }
  library("Matrix")

  if("doParallel" %in% rownames(installed.packages()) == FALSE)
  {install.packages("doParallel")
  }
  #library("doParallel")

  if("parallel" %in% rownames(installed.packages()) == FALSE)
  {install.packages("parallel")
  }
  library("parallel")

  if("splines" %in% rownames(installed.packages()) == FALSE)
  {install.packages("splines")
  }
  library("splines")
}
