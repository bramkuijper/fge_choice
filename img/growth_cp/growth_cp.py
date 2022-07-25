#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import re


# function to find the start of the parameter
# part in the data file
def find_parameter_line(filename):
    f = open(filename)
    fl = f.readlines();
    f.close()

    # make a reverse sequence (shorter looping)
    revseq = range(len(fl) - 1,0,-1)

    # search for first line that starts with a number
    # and finish there
    for line_i in revseq:
        if re.search("^\d.*",fl[line_i]) != None:
            return(line_i)

    return(None)

# set basic theme
sns.set_theme(style="ticks")

# file name 
file_choosy_outperforms = "choosy_outperforms.csv"

# get line at which data file ends
line_data_ends = find_parameter_line(
        filename = file_choosy_outperforms)

# the data in which choosy outperforms anything else
data_choosy_outperforms = pd.read_csv(
        filepath_or_buffer=file_choosy_outperforms
        ,sep=";"
        ,nrows=line_data_ends - 1)

data_choosy_outperforms["Promiscuous"] = data_choosy_outperforms["Ipg1"] + data_choosy_outperforms["Ipg2"]
data_choosy_outperforms["Choosy"] = data_choosy_outperforms["Icg2"] + data_choosy_outperforms["Icg1"]

# select columns we would like to plot
subset_choosy = data_choosy_outperforms[["time","Choosy","Promiscuous"]]

# data currently in wide format, transform to long
subset_choosy_long = pd.melt(frame=subset_choosy 
        ,id_vars=["time"]
        ,value_vars=["Choosy","Promiscuous"])


the_fig = sns.relplot(
        data=subset_choosy_long.loc[subset_choosy_long["time"] < 20000]
        ,x="time"
        ,y="value"
        ,hue="variable"
        ,kind="line"
        )

the_fig

the_fig.savefig(fname="plot_choosy_vs_promiscuous.pdf")
