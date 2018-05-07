load_data <- function() {
  list_files <- list.files('../data', full.names = TRUE, recursive = TRUE )
  list_files <- list_files[grepl('.Rdata', list_files, ignore.case = FALSE)]
  
  mylist<- lapply(list_files, function(x) {
    load(file = x)
    get(ls()[ls()!= "filename"])
  })
  
  names(mylist) <- list_files
  mylist
  #list2env(mylist ,.GlobalEnv)
  }




