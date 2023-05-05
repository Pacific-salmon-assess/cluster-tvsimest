#' get slurm results out without a aslurm job object
my_get_slurm_out <- function (slr_job_name, nodes.list=0:50, outtype = "raw") 
{
       
    res_files <- paste0("results_", nodes.list, ".RDS")
    tmpdir <- paste0("_rslurm_", slr_job_name)
    missing_files <- setdiff(res_files, dir(path = tmpdir))
    if (length(missing_files) > 0) {
        missing_list <- paste(missing_files, collapse = ", ")
        warning(paste("The following files are missing:", missing_list))
    }
    res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
    if (length(res_files) == 0) 
        return(NA)
    
        slurm_out <- lapply(res_files, readRDS)
    
    slurm_out <- do.call(c, slurm_out)
    if (outtype == "table") {
        slurm_out <- as.data.frame(do.call(rbind, slurm_out))
    }
    return(slurm_out)
}

