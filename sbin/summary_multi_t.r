#!/usr/bin/env Rscript 

library("tidyverse")
library("here")
library("parallel")


# obtain the final line of the actual simulation data
# of each file
get_last_line_number <- function(file_name)
{
    # get all the lines
    all_lines <- readLines(con=file_name)
    
    # obtain number of lines
    nlines <- length(all_lines)
    
    # search for last line from the end
    # (coz faster)
    for (line_idx in seq(nlines,1,-1))
    {
        if (length(grep(pattern="^\\d",all_lines[[line_idx]]))>0)
        {
            return(line_idx)
        }
    }
    
    stop("no line found")
}

# function to obtain the parameter listing at the bottom
# of each file
get_parameters <- function(file_name)
{
    last_line <- get_last_line_number(file_name = file_name)
    
    f <- read.table(file = file_name
                    ,skip = last_line
                    ,sep=";"
                    ,header=F)
    
    # get row with values 
    vals <- as.data.frame(t(f[,2]))
    
    # add names 
    names(vals) <- f[,1]
    
    return(vals)
}

extract_from_single_file <- function(file_name, timepoints)
{
    # get the last line of this file
    last_line_i = get_last_line_number(file_name = file_name)

    data <- read_delim(
            file=file_name
            ,n_max=last_line_i - 1
            ,show_col_types=F
            ,delim=";"
            )

    # obtain final timepoint (we would also like to include
    # it to check what happens at equilibrium)
    timepoint_last <- data[nrow(data),"time"]

    all_time_points <- c(timepoints, timepoint_last)

    # get data at desired time point
    data_timepoints <- filter(data, time %in% all_time_points)
    data_timepoints[,"file"] <- file_name

    return(data_timepoints)
} # end extract_from_single_file

# get the data at time point tdata
get_data <- function(path, tdata, filename_regexp)
{
    all_files <- list.files(
            path=path
            ,pattern=filename_regexp
            ,recursive = F
            ,full.names=T)

    stopifnot(length(all_files) > 0)

    all_data <- do.call(
            rbind,
            mclapply(X = all_files
            ,FUN = extract_from_single_file
            ,timepoints=tdata
            ,mc.cores=20L)
    )

    return(all_data)
} # end get_data()


path = commandArgs(trailingOnly=T)[[1]]
pattern = commandArgs(trailingOnly=T)[[2]]

# get the data from the numeric solver
data_numeric_solver <- get_data(
        path=path,
        tdata=c(0,1,5,10,50,100,1000),
        filename_regexp=pattern
        )

# make delta data
data_numeric_solver <- mutate(
        data_numeric_solver
        ,pBc=ICG1/(SC + ICG1 + ICG2)
        ,pBp=IPG1/(SP + IPG1 + IPG2)
        ,pGc=ICG2/(SC + ICG1 + ICG2)
        ,pGp=IPG2/(SP + IPG1 + IPG2)
        ,fc=(ICG1 + ICG2 + SC)/N 
        )

write_delim(x=data_numeric_solver
        ,delim=";"
       ,file="summary.csv")
