library("tidyverse")
library("here")
library("patchwork")
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
    } else if (dataset[1,"IPG1"] < 1 && dataset[1,"IPG2"] > 1)
    {
        single_or_mixed <- "M13d"
        order <- 3
    } else if (dataset[1,"IPG1"] > 1 && dataset[1,"IPG2"] < 1)
    {
        single_or_mixed <- "M13s"
        order <- 2
    }
    
    dataset[,"single_or_mixed"] <- single_or_mixed
    dataset[,"order"] <- order

    return(dataset)
} # end find_out_type()

file_names <- c("output2","output3","output4")

label <- c("M13s_only","M13d_only","Mixed")


for (file_idx in 1:length(file_names))
{
    file_time_course_antibiotics = paste0(
            "../selection/base_output_low_advantage_M13s/", file_names[[file_idx]]
    )

    last_line_i = get_last_line_number(file_name = file_time_course)
        
    time_series_data <- read_delim(file = path
                                   ,delim =";"
                                   ,n_max = last_line_i - 1)


    clip <- function(x,xmin=0,xmax=1)
    {
        return(ifelse(x < xmin, xmin, ifelse(x > xmax, xmax, x)))
    }



    time_series_data <- time_series_data %>% mutate(
        M13s = IPG1 + ICG1
        ,M13d = IPG2 + ICG2
        ,fc = clip((ICG1 + ICG2 + SC) / (ICG1 + ICG2 + SC + IPG1 + IPG2 + SP),xmin=0.0,xmax=1.0)
    )

    time_series_data_first <- time_series_data %>% filter(time < 1000)


    # transform so that we have M13d and M13s in one plot
    time_series_data_l <- pivot_longer(data=time_series_data_first
                                       ,cols=c(M13s,M13d)
                                       ,values_to = "Frequency"
                                       ,names_to = "FGE_type")

    # what do we want to plot: M13d M13s
    fge_time <- ggplot(data=time_series_data_l
           ,mapping=aes(x=time
                        ,y=Frequency)) +
       geom_line(mapping=aes(colour=FGE_type),linewidth=1) +
        theme_classic(base_size = 18) +
        scale_colour_brewer(palette = "Set1") +
        xlab("Time step")




    # what do we want to plot: M13d M13s
    fc_time <- ggplot(data=time_series_data_l
           ,mapping=aes(x=time
                        ,y=fc)) +
       geom_line(linewidth=1) +
        theme_classic(base_size = 18) +
        scale_colour_brewer(palette = "Set1") +
        ylim(0,1) + 
        xlab("Time step") +
        ylab("Frequency of CI allele")

    fge_time | fc_time

    ggsave(paste0("time_course_FGE",label[[file_idx]],".pdf")
            ,width=10)

} # end for file_idx
