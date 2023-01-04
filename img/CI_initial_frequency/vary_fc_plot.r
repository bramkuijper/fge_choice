library("tidyverse")
library("here")
library("cowplot")
library("patchwork")

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
    for (file_i in all_files)
    {
        last_line_i = get_last_line_number(file_name = file_i)
    
        data <- read_delim(file=file_i, n_max=last_line_i - 1)


        # get data at desired time point
        data_timepoints <- filter(data, time %in% tdata)
        data_timepoints[,"file"] <- file_i

        all_data <- bind_rows(all_data, data_timepoints)
    }

    return(all_data)
} # end get_data()

########## PATHS, PARAMETERS, ETC ##########

main_path = here()

path=file.path(main_path, "img/CI_initial_frequency/data_original")




if (!exists("data_numeric_solver"))
{
    # get the data from the numeric solver
    data_numeric_solver <- get_data(
            path=path,
            tdata=c(0,1,5,10,50,100,1000),
            filename_regexp="^varyfc*"
            )
}

# make delta data
data_numeric_solver <- mutate(
        data_numeric_solver
        ,pBc=ICG1/(SC + ICG1 + ICG2)
        ,pBp=IPG1/(SP + IPG1 + IPG2)
        ,pGc=ICG2/(SC + ICG1 + ICG2)
        ,pGp=IPG2/(SP + IPG1 + IPG2)
        ,fc=(ICG1 + ICG2 + SC)/N 
        ,fM13s=(ICG1 + IPG1)/N 
        ,fM13d=(ICG2 + IPG2)/N 
        )

# then calculate differences for each series of times
# see https://stackoverflow.com/questions/72534996/r-subtract-rows-from-the-first-row-in-each-group 
data_numeric_solver_t <- data_numeric_solver %>% group_by(file) %>% arrange(time) %>%
        mutate(
            delta_pGc = pGc - pGc[row_number() == 1]
            ,delta_pGp = pGp - pGp[row_number() == 1]
            ,delta_pBc = pBc - pBc[row_number() == 1]
            ,delta_pBp = pBp - pBp[row_number() == 1]
            ,delta_M13s = fM13s - fM13s[row_number() == 1]
            ,delta_M13d = fM13d - fM13d[row_number() == 1]
            ,init_fc=fc[row_number() == 1]
            )

# make two tibbles in long format: one with the change in conditional 
# frequencies, the other with the change in overall freq of M13s vs M13d
data_numeric_solver_tl <- data_numeric_solver_t %>% filter(time > 0) %>%
    pivot_longer(
        ,cols=c(delta_pGc, delta_pBc, delta_pGp, delta_pBp)
        ,names_to="conditional_type"
        ,values_to="conditional_value") 

data_numeric_solver_tl_M13 <- data_numeric_solver_t %>% filter(time > 0) %>%
    pivot_longer(
        ,cols=c(delta_M13s, delta_M13d)
        ,names_to="FGE_type"
        ,values_to="FGE_freq") 


# first global plot of conditional change over x timesteps
ggplot(data=data_numeric_solver_tl
        ,mapping=aes(x=init_fc, y = conditional_value)) +
    geom_line(mapping=aes(colour=conditional_type)) +
    facet_wrap(vars(time), ncol=4) +
    scale_colour_manual(name=""
            ,values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
            ,breaks=c("delta_pGc","delta_pBc","delta_pGp","delta_pBp")
            ,labels=c("M13d | CI", "M13s | CI", "M13d | S", "M13s | S")) +
    theme_bw() +
    xlab("Initial frequency of CI") +
    ylab("Change in frequency of M13") +
    ylim(-0.1,1)

ggsave("vary_fc_across_time_points.pdf")


# patchwork plot of the two figs 
# spread out over CI or total 
p1 <- ggplot(data=filter(data_numeric_solver_tl_M13, time ==  100)
        ,mapping=aes(x=init_fc, y = FGE_freq)) + 
    geom_line(mapping=aes(colour=FGE_type),size=1) +
    theme_classic(base_size=22) +
    xlab("Initial frequency of CI") +
    ylab("Change in frequency of M13") +
    ylim(-0.05,1.0)

p2 <- ggplot(data=filter(data_numeric_solver_tl,time == 100)
        ,mapping=aes(x=init_fc, y = conditional_value)) +
    geom_line(mapping=aes(colour=conditional_type),size=1) +
    scale_colour_manual(name=""
            ,values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
            ,breaks=c("delta_pGc","delta_pBc","delta_pGp","delta_pBp")
            ,labels=c("M13d | CI", "M13s | CI", "M13d | S", "M13s | S")) +
    theme_classic(base_size=22) +
    xlab("Initial frequency of CI") +
    ylab("Change in frequency of M13") +
    ylim(-0.05,1)


p1 | p2

ggsave("vary_fc_timepoint_100.pdf",width=15,height=5)





