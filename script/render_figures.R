for (n in c(1,2,3,5,6)) {
  rmarkdown::render(sprintf("fig%d/fig%d.Rmd", n,n), envir = new.env(parent = globalenv()))
}
