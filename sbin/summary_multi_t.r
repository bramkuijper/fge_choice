#!/usr/bin/env Rscript 

library("tidyverse")
library("here")


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

# get the data at time point tdata
get_data <- function(path, tdata, filename_regexp)
{
    all_files <- list.files(
            path=path
            ,pattern=filename_regexp
            ,recursive = F
            ,full.names=T)

    stopifnot(length(all_files) > 0)

    # allocate empty variable
    all_data <- NULL

    # get last line
    for (file_idx in 1:length(all_files))
    {
        file_i <- all_files[[file_idx]]
    
        print(paste0("File ",file_idx," out of ",length(all_files)))

        last_line_i = get_last_line_number(file_name = file_i)
    
        data <- read_delim(
                file=file_i 
                ,n_max=last_line_i - 1
                ,show_col_types=F
                )


        # get data at desired time point
        data_timepoints <- filter(data, time %in% tdata)
        data_timepoints[,"file"] <- file_i

        all_data <- bind_rows(all_data, data_timepoints)
    }

    return(all_data)
} # end get_data()


path = commandArgs(trailingOnly=T)[[1]]

# get the data from the numeric solver
data_numeric_solver <- get_data(
        path=path,
        tdata=c(0,1,5,10,50,100,1000),
        filename_regexp="^varyfc*"
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

