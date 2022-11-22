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


##### PLOT WITH ANTIBIOTICS #####

path=file.path(main_path, "/img/selection/base_output_low_advantage_M13s/")


# get the data from the numeric solver
data_numeric_solver <- get_data(
        path=path, 
        tdata=c(0,50,100, 500, 1000, 5000, 9000),
        filename_regexp="^output_.*"
)

# calculate conditional frequencies
data_numeric_solver <- mutate(data_numeric_solver,
        # p_{B \mid c} = I_{cB} / (S_{c} + I_{cB} + I_{cG})
#        pBc=ICG1/(SC + ICG1 + ICG2)
#        ,pBp=IPG1/(SP + IPG1 + IPG2)
#        ,pGc=ICG2/(SC + ICG1 + ICG2)
#        ,pGp=IPG2/(SP + IPG1 + IPG2)
        pTotal=(IPG1 + IPG2 + SP) / N
        ,cTotal=(ICG1 + ICG2 + SC) / N
        )

# calculate deltas
data_numeric_solver_delta = 
    data_numeric_solver %>% group_by(file) %>% arrange(time) %>% mutate(
        delta_p = pTotal - pTotal[row_number() == 1]
        ,delta_c = cTotal - cTotal[row_number() == 1]
    )

ggplot(data=data_numeric_solver_delta %>% filter(time==500)
        ,mapping=aes(x=fct_reorder(single_or_mixed,order)
                     ,y=delta_c)) + 
    geom_bar(stat="identity") +
    theme_classic(base_size = 18) +
    xlab("") +
    ylab("Selection rate")

ggsave(file="selection_with_antibiotics.pdf")

#################### ANTIBIOTICS LOW COST OF CHOICE ####################

rm(data_numeric_solver)

path=file.path(main_path, "img/selection/low_c_output/")


# get the data from the numeric solver
data_numeric_solver <- get_data(
        path=path, 
        tdata=c(0,50,100, 500, 1000, 5000, 9000),
        filename_regexp="^output_.*"
)

# calculate conditional frequencies
data_numeric_solver <- mutate(data_numeric_solver,
        # p_{B \mid c} = I_{cB} / (S_{c} + I_{cB} + I_{cG})
#        pBc=ICG1/(SC + ICG1 + ICG2)
#        ,pBp=IPG1/(SP + IPG1 + IPG2)
#        ,pGc=ICG2/(SC + ICG1 + ICG2)
#        ,pGp=IPG2/(SP + IPG1 + IPG2)
        pTotal=(IPG1 + IPG2 + SP) / N
        ,cTotal=(ICG1 + ICG2 + SC) / N
        )

# calculate deltas
data_numeric_solver_delta = 
    data_numeric_solver %>% group_by(file) %>% arrange(time) %>% mutate(
        delta_p = pTotal - pTotal[row_number() == 1]
        ,delta_c = cTotal - cTotal[row_number() == 1]
    )

ggplot(data=data_numeric_solver_delta %>% filter(time==500)
        ,mapping=aes(x=fct_reorder(single_or_mixed,order)
                     ,y=delta_c)) + 
    geom_bar(stat="identity") +
    theme_classic(base_size = 18) +
    xlab("") +
    ylab("Selection rate")

ggsave(file="selection_with_antibiotics_low_cost.pdf")

#################### NO ANTIBIOTICS ####################

rm(data_numeric_solver)

path=file.path(main_path, "img/selection/no_antibiotics/")


# get the data from the numeric solver
data_numeric_solver <- get_data(
        path=path, 
        tdata=c(0,50,100, 500, 1000, 5000, 9000),
        filename_regexp="^output_.*"
)

# calculate conditional frequencies
data_numeric_solver <- mutate(data_numeric_solver,
        # p_{B \mid c} = I_{cB} / (S_{c} + I_{cB} + I_{cG})
#        pBc=ICG1/(SC + ICG1 + ICG2)
#        ,pBp=IPG1/(SP + IPG1 + IPG2)
#        ,pGc=ICG2/(SC + ICG1 + ICG2)
#        ,pGp=IPG2/(SP + IPG1 + IPG2)
        pTotal=(IPG1 + IPG2 + SP) / N
        ,cTotal=(ICG1 + ICG2 + SC) / N
        )

# calculate deltas
data_numeric_solver_delta = 
    data_numeric_solver %>% group_by(file) %>% arrange(time) %>% mutate(
        delta_p = pTotal - pTotal[row_number() == 1]
        ,delta_c = cTotal - cTotal[row_number() == 1]
    )
print(data_numeric_solver_delta %>% filter(time==500))
ggplot(data=data_numeric_solver_delta %>% filter(time==500)
        ,mapping=aes(x=fct_reorder(single_or_mixed,order)
                     ,y=delta_c)) + 
    geom_bar(stat="identity") +
    theme_classic(base_size = 18) +
    xlab("") +
    ylab("Selection rate") 

ggsave(file="selection_no_antibiotics.pdf")
