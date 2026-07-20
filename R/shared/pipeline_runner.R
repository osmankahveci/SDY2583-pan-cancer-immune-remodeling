# Ordered panel-pipeline execution helper.
run_sdy2583_script_sequence <- function(panel, scripts, validate_script=NULL, stop_on_missing=TRUE) {
  root<-normalizePath(path.expand(Sys.getenv("SDY2583_REPO_ROOT",unset=getwd())),mustWork=FALSE)
  Sys.setenv(SDY2583_REPO_ROOT=root)
  config<-file.path(root,"config","paths.R"); if(file.exists(config))source(config)
  status<-data.frame(panel=character(),script=character(),status=character(),message=character(),stringsAsFactors=FALSE)
  for(relative in scripts){
    path<-file.path(root,relative)
    if(!file.exists(path)){
      status<-rbind(status,data.frame(panel=panel,script=relative,status="missing",message="file not found"))
      if(stop_on_missing)stop("Missing pipeline script: ",path) else next
    }
    cat("\n============================================================\n",panel," :: ",relative,"\n============================================================\n",sep="")
    result<-tryCatch({source(path,local=new.env(parent=globalenv()));list(ok=TRUE,message="completed")},error=function(e)list(ok=FALSE,message=conditionMessage(e)))
    status<-rbind(status,data.frame(panel=panel,script=relative,status=ifelse(result$ok,"completed","failed"),message=result$message))
    if(!result$ok)stop(panel," pipeline failed in ",relative,": ",result$message)
  }
  if(!is.null(validate_script)){
    path<-file.path(root,validate_script)
    if(file.exists(path))source(path,local=new.env(parent=globalenv())) else if(stop_on_missing)stop("Validation script missing: ",path)
  }
  invisible(status)
}
