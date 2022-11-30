library("tidyverse")
path <- "."

all_files <- list.files(
        path=path
        ,pattern="^output_.*"
        ,recursive = F
        ,full.names=T)

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

# get all the parameters across all files
if (!exists("all_data"))
{
    all_data <- NULL
    
    stopifnot(length(all_files) > 0)
    
    # go through all files and get parameters
    # and collect them in one data frame
    for (file in all_files)
    {
        # get the last line of the data set
        last_line_data <- get_last_line_number(file_name = file)
        
        current_data <- read_delim(file=file,delim=";",n_max = last_line_data - 1)
        
        # add file name
        current_data[,"type"] <- file

        # make names uppercase
        names(current_data) <- toupper(names(current_data))


        
        label <- "Both MGEs"
        order <- "d"
        
        # find out what MGE we are talking about here (this depends
        # on initial conditions, note 1st row subset select) and label accordingly
        # here Ipg1 is promiscuous bad FGE, Ipg2 is promiscuous good FGE
        # Icg1 is choosy bad FGE, Ipg2 is choosy good FGE
        if (current_data[1,"IPG1"] < 1 && current_data[1,"IPG2"] < 1)
        {
            label <- "No MGE"
            order <- "a"
        } else if (current_data[1,"IPG1"] < 1 && current_data[1,"IPG2"] > 1)
        {
            label <- "Only MGE G"
            order <- "c"
        } else if (current_data[1,"IPG1"] > 1 && current_data[1,"IPG2"] < 1)
        {
            label <- "Only MGE B"
            order <- "b"
        }
        
        current_data[,"label"] <- label
        current_data[,"order"] <- order
        
        # concat data
        all_data <- bind_rows(all_data,current_data)
    }
}


tmax <- 20000

# subset data dependent on tmax and then pivot_longer
# to make sure we can plot lines
all_data <- all_data %>% filter(TIME < tmax) %>% mutate(
    choosy=ICG1 + ICG2 
    ,promiscuous=IPG1+IPG2
    ,time_thousands = TIME/1000
)

all_data_t <- all_data %>% pivot_longer(
    cols=c("choosy","promiscuous")
    ,names_to="Type"
    ,values_to="Frequency"
)

# some labelling vector
# really not sure why it has to be a named
# character vector, simply order of the panels
# should be enough here, but Hadley says no
the_labels <- c(
    a="No MGE"
    ,b="Only MGE B"
    ,c="Only MGE G"
    ,d="Both MGEs"
)

ggplot(data=all_data_t
       ,mapping=aes(x=time_thousands
                    ,y=Frequency)) +
    geom_line(mapping=aes(colour=Type)) +
    scale_fill_brewer(palette = "Set2") +
    theme_classic() +
    facet_grid(~order, labeller = labeller(order=the_labels)) +
    scale_colour_manual(values=c("#92d4f2","#992f67"), labels=c("Immune","Sensitive")) +
    xlab("Time step (x1000)") +
    ylab("Density")

ggsave(file="choosy_promiscuous_vs_time.pdf"
       ,height=3)


