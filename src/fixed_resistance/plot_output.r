#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}

xvar <- "time"

jsonstuff <- paste0('[
    {"xvar" : "',xvar,'",
    "yvar" : ["Sp","Sc"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["Ipg1","Ipg2","Icg1","Icg2"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : "N"
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["Choosy","Promiscuous"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["log10Choosy","log10Promiscuous"]
    }
]')

if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

data.tibble.orig <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

data.tibble.orig <- mutate(data.tibble.orig,
        Ipg2 = Ipg2 
        ,Choosy=Icg1 + Icg2
        ,Promiscuous=Ipg1 + Ipg2
        ,log10Choosy=log10(Icg1 + Icg2)
        ,log10Promiscuous=log10(Ipg1 + Ipg2)
        )

if (nrow(data.tibble.orig) > 50000)
{
    data.tibble <- filter(data.tibble.orig, time %% 100 == 0)
} else
{
    data.tibble <- data.tibble.orig
}



plot.structure <- fromJSON(jsonstuff, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list(rep(NA,times=plot.structure.l*2))

plot.list.idx <- 1

for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes_string(x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values)) + geom_line(mapping=aes_string(colour=yvar_name))
    } else {
        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble
                ,mapping=aes_string(x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list.idx <- plot.list.idx + 1
}

data.tibble.sub <- data.tibble.orig[data.tibble.orig$time < 2000,]

for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble.sub
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes_string(x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values)) + geom_line(mapping=aes_string(colour=yvar_name))
    } else {
        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble.sub
                ,mapping=aes_string(x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list.idx <- plot.list.idx + 1
}

wrap_plots(plot.list,ncol=1)

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height= 3 * plot.structure.l)
