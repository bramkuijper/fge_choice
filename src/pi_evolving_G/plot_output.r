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

jsonstuff <- '[
    {"xvar" : "time_step",
    "yvar" : ["S","G1","G2","G1G2"]
    },
    {
        "xvar" : "time_step",
        "yvar" : ["uS","uG1","uG2","uG1G2"]
    },
    {
        "xvar" : "time_step",
        "yvar" : ["vS","vG1","vG2","vG1G2"]
    },
    {
        "xvar" : "time_step",
        "yvar" : "N"
    },
    {
        "xvar" : "time_step",
        "yvar" : "pi"
    },
    {
        "xvar" : "time_step",
        "yvar" : "delta_pi"
    }
]
'

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

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

if (nrow(data.tibble) > 50000)
{
    data.tibble <- filter(data.tibble, time_step %% 100 == 0)
}

plot.structure <- fromJSON(jsonstuff, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list(rep(NA,times=plot.structure.l))

for (plot_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        plot.list[[plot_idx]] <- ggplot(data=sub.data
                ,mapping=aes_string(x=plot.structure[[plot_idx]]$xvar
                        ,y=yvar_values)) + geom_line(mapping=aes_string(colour=yvar_name))
    } else {
        plot.list[[plot_idx]] <- ggplot(data=data.tibble
                ,mapping=aes_string(x=plot.structure[[plot_idx]]$xvar
                        ,y=plot.structure[[plot_idx]]$yvar)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_idx]]))
    {
        plot.list[[plot_idx]] <- plot.list[[plot_idx]] + ylim(
                unlist(
                        plot.structure[[plot_idx]]$ylim)
                )
    }
}

wrap_plots(plot.list,ncol=1)

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height= 3 * plot.structure.l)






