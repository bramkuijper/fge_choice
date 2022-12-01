# ratio plot, showing the ratio
# of infection between immune and sensitive 

library("tidyverse")
library("here")
library("cowplot")
library("forcats")

main_path = here()


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

# find out what the type is of this data
find_out_type <- function(dataset)
{
    single_or_mixed <- "mixed"
    
    order <- 4

    if (dataset[1,"IPG1"] < 1 && dataset[1,"IPG2"] < 1)
    {
        single_or_mixed <- "none"
        order <- 1
    } else if (dataset[1,"IPG1"] < 1 && dataset[1,"IPG2"] >= 1)
    {
        single_or_mixed <- "M13d"
        order <- 3
    } else if (dataset[1,"IPG1"] >= 1 && dataset[1,"IPG2"] < 1)
    {
        single_or_mixed <- "M13s"
        order <- 2
    }
    
    dataset[,"single_or_mixed"] <- single_or_mixed
    dataset[,"order"] <- order

    return(dataset)
} # end find_out_type()


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
    for (file_i in all_files)
    {
        last_line_i = get_last_line_number(file_name = file_i)
    
        data <- read_delim(file=file_i, n_max=last_line_i - 1)

        # determine type
        data <- find_out_type(data)    

        data[,"file"] <- file_i

        all_data <- bind_rows(all_data, data[data$time %in% tdata,])
    }

    return(all_data)
} # end get_data()

path = "../time_course_facets/"

# first get the summary file
summary_data <- read_delim(
        file=file.path(path,"summary.csv")
        ,delim=";")

# find an antibiotic file that matches the proper value of psi and FG1
file_name_antibiotic = summary_data %>% filter(psiG1 == 1 & FG1 == 2 & FG2 == 10) %>% pull(file)


file_name_no_antibiotic = summary_data %>% filter(psiG1 == 1 & FG2 == 1) %>% pull(file)

t_measure <- 750

data_numeric_solver_no_antibiotic <- read_delim(
        file=file.path(path,file_name_no_antibiotic)
        ,delim=";"
        ,n_max = t_measure + 1) %>% filter(time == t_measure)


data_numeric_solver <- read_delim(
        file=file.path(path,file_name_antibiotic)
        ,delim=";"
        ,n_max = t_measure + 1) %>% filter(time == t_measure)

# we need to plot ratios of mixed infections 
path=file.path(main_path, "/img/time_course_facets/")



# infection ratio immune vs sensitive M13s
M13s_ratio_no_anti <- data_numeric_solver_no_antibiotic$ICG1 / data_numeric_solver_no_antibiotic$IPG1
M13d_ratio_no_anti <- data_numeric_solver_no_antibiotic$ICG2 / data_numeric_solver_no_antibiotic$IPG2
M13s_ratio <- data_numeric_solver$ICG1 / data_numeric_solver$IPG1
M13d_ratio <- data_numeric_solver$ICG2 / data_numeric_solver$IPG2

# now build a box plot of the ratios
the.ratios <- data.frame(
        x=c(1,2,3,4)
        ,y=c(M13s_ratio_no_anti
                ,M13d_ratio_no_anti
                ,M13s_ratio
                ,M13d_ratio)
        ,labels=c("M13s","M13d","M13s","M13d")
    )


ggplot(data=the.ratios
        ,mapping=aes(x=factor(x)
                ,y=y)) +
    geom_point(size=5,colour="#337167") +
    theme_classic(base_size=18) +
    xlab("") +
    ylab("Ratio of infection\nimmune/sensitive") +
    geom_hline(yintercept=1.0,linetype="dashed") +
    scale_x_discrete(breaks=waiver(), labels=the.ratios$labels)  +
    scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1.0,1.25)) 

ggsave("ratio_plot.pdf")

