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

single_multipanel <- function(
        file_antibiotic
        ,file_no_antibiotic
        ,file_out
        ,tmax=200)
{
    file_names <- c(file_antibiotic,
                    file_no_antibiotic)
    
    antibiotics <- c("antibiotics","no_antibiotics")
    
    all_time_series_data <- NULL
    
    path ="."
    
    # collect the data from the files
    for (file_idx in 1:length(file_names))
    {
        file_time_course = file.path(
                path, file_names[[file_idx]]
        )
        
       last_line_i = get_last_line_number(file_name = file_time_course)
            
        time_series_data <- read_delim(file = file_time_course
                                       ,delim =";"
                                       ,n_max = last_line_i - 1)
    
    
        clip <- function(x,xmin=0,xmax=1)
        {
            return(ifelse(x < xmin, xmin, ifelse(x > xmax, xmax, x)))
        }
    
        time_series_data <- time_series_data %>% mutate(
            M13s = IPG1 + ICG1
            ,M13d = IPG2 + ICG2
            ,ImmuneM13s= (ICG1) / (ICG1 + ICG2 + SC)
            ,ImmuneM13d= (ICG2) / (ICG1 + ICG2 + SC)
            ,SensitiveM13s= (IPG1) / (IPG1 + IPG2 + SP)
            ,SensitiveM13d= (IPG2) / (IPG1 + IPG2 + SP)
            ,fc = clip((ICG1 + ICG2 + SC) / (ICG1 + ICG2 + SC + IPG1 + IPG2 + SP),xmin=0.0,xmax=1.0)
            ,antibiotics= antibiotics[[file_idx]]
        )
    
        time_series_data_first <- time_series_data %>% filter(time < tmax)
    
        all_time_series_data <- bind_rows(
                all_time_series_data, 
                time_series_data_first)
    
    } # end for
    
    parameters <- get_parameters(file_name = file_antibiotic)
    
    # we need to get a single column for all the frequency data
    all_time_series_data_l = all_time_series_data %>% pivot_longer(
        cols = c("ImmuneM13s","ImmuneM13d","SensitiveM13s","SensitiveM13d")
        ,values_to = "Frequency"
        ,names_to = "Type"
    )
    
    # then we need to split the frequency data up according to a panel
    all_time_series_data_l <- all_time_series_data_l %>% mutate(
                m13=ifelse(regexpr(pattern="M13s",text=Type)!=-1,
                           "M13s","M13d")
                ,ctype=ifelse(regexpr(pattern="Immune",text=Type)!=-1,
                              "Immune","Susceptible")
                )
    
    # now make the plot
    ggplot(data=all_time_series_data_l
            ,mapping=aes(x=time
                         ,y=Frequency)) +
        geom_line(mapping = aes(colour=ctype)) +
    #    facet_grid(~ fct_relevel(label,"M13s_only","M13d_only")) +
        facet_grid(antibiotics~m13) +
        theme_classic(base_size=16) +
        scale_colour_brewer(palette = "Set1") +
        xlab("Time") +
        ylab("Frequency") +
        ylim(0,1) + plot_annotation(
        title=paste0("FG1: "
            ,parameters["FG1"]))
    
    ggsave(file=file_out) 
} # end single_multipanel


# first read in the summary of all the data files
summary_data <- read_delim(file = "summary_iterations.csv"
                           ,delim = ";")

# first get the antibiotic data
summary_data_antibiotic <- summary_data %>% filter(FG1 > 1)

summary_data_no_antibiotic <- summary_data %>% filter(FG1 == 1)

# should be a single line of no antibiotic data
stopifnot(nrow(summary_data_no_antibiotic) == 1)

for (row_i in 1:nrow(summary_data_antibiotic))
{
    file_antibiotics <- summary_data_antibiotic[row_i,] %>% pull(file)
    
    stopifnot(nrow(file_antibiotics) == 1)
    
    single_multipanel(
            file_antibiotic=file_antibiotics
            ,file_no_antibiotic=summary_data_no_antibiotic %>% pull(file)
            ,file_out=paste0(basename(file_antibiotics),".pdf")
            ,tmax=1000)
}