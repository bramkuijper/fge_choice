library("tidyverse")
library("here")
library("cowplot")

main_path = here()

path=file.path(main_path, "src/fixed_resistance/")

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

get_data <- function(path)
{
    all_files <- list.files(
            path=path
            ,pattern="^output_.*"
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

        data[,"file"] <- file_i

        all_data <- bind_rows(all_data, data[nrow(data),])
    }

    return(all_data)
} # end get_data()

the_data <- get_data(path=path)

# calculate conditional frequencies
the_data <- mutate(the_data,


# make first panel: p and g vs frequency of infection with M13

bar_padding <- 0.1
bar_width <- 0.5

# simulate some data first to get the plotting correct
group_vals <- c(0,0,1,1,2.25,2.25,3.25,3.25)

the_data=data.frame(
        group=group_vals,
        host_type=factor(
                rep(c("S","CI"),length.out=8)))

# colour values
colors <- c("#88ccee","#882255")


# make list with x label positions. Totally no clue how gglplot calculates the 
# actual distances between padded bars. This is all a shot in the dark
x_label_pos_delta <- 0.25 * rep(c(-(bar_width-bar_padding),(bar_width-bar_padding)),length.out=length(group_vals))
x_label_pos_abs <- group_vals + x_label_pos_delta

# make list with x labels
x_labels <- rep(c("S","CI"),length.out=length(group_vals))

g1 <- ggplot(data=the_data
        ,mapping=aes(x=group
                ,y=yval,
                ,fill=host_type)) +
    geom_bar(stat="identity"
            ,width=bar_width
            ,position=position_dodge2(
                    preserve="single"
                    ,padding=padding_bar)
            ,colour="black") +
    geom_vline(mapping=aes(xintercept=1+(2.25-1.0)/2.0)
            ,colour="black"
            ,linetype=2
            ) +
    geom_vline(mapping=aes(xintercept=0.5)
            ,colour="grey"
            ,linetype=3
            ) +
    geom_vline(mapping=aes(xintercept=2.75)
            ,colour="grey"
            ,linetype=3
            ) +
    scale_x_continuous(name=""
            ,breaks=x_label_pos_abs
            ,labels=x_labels
            ,limits=c(min(x_label_pos_abs)-0.75,max(x_label_pos_abs)+0.75)
            ,expand=expansion(mult=c(0.1,0.1))) +

    #    scale_y_continuous(name=""
    #            ,limits=c(0,1)
    #            ,expand=expansion(mult=c(0.01,0.01))) +

    # manual colour scheme
    scale_fill_manual(values = colors) +
    coord_cartesian(ylim = c(0, 1)
            ,clip="off"
            ,xlim=c(min(group_vals),max(group_vals))) +
    theme_classic() +

    # some mods on theme: black labels + no legend
    theme(axis.text.x = element_text(color="black") # make tick label text black
            ,axis.text.y = element_text(color="black")
            ,legend.position="none") +

    # infection label, left bottom
    annotate(geom="text"
            ,x=c(-0.7,sort(unique(group_vals)))
            ,y=-0.1
            ,label=c("infection:",rep(c("M13s","M13d"),times=2))) +
    
    # mixed infection label, top right
    annotate(geom="text"
            ,x=group_vals[[6]] + (group_vals[[8]] - group_vals[[6]])/2.0
            ,y=0.9
            ,label=c("Mixed M13"))

gt <- ggplot_gtable(ggplot_build(g1))
gt$layout$clip[gt$layout$name=="panel"] <- "off"

ggsave("freq_w_M13.pdf" , gt)


