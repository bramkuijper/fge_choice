#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import re

from matplotlib import rcParams
import matplotlib as mpl

fontpath = "/System/Library/Fonts/Supplemental/" 
fontpath = "/usr/local/texlive/2022/texmf-dist/fonts/opentype/public/stix2-otf/"

mpl.use("pgf")

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "pgf.texsystem": "xelatex",
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
        r"\usepackage{units}",         # load additional packages
        r"\usepackage{mathspec}",         # load additional packages
        r"\setmainfont[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = *-Regular," +\
            "ItalicFont = *-Italic," +\
            "BoldFont = *-Bold," +\
            "Extension = .otf]{STIXTwoText}",
        r"\setmathfont[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = *-Regular," +\
            "ItalicFont = *-Italic," +\
            "BoldFont = *-Bold," +\
            "Extension = .otf]{STIXTwoText}",
        r"\setmathrm[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = *-Regular," +\
            "ItalicFont = *-Italic," +\
            "BoldFont = *-Bold," +\
            "Extension = .otf]{STIXTwoText}",
         ])
}

mpl.rcParams.update(pgf_with_custom_preamble)


#pgf_with_custom_preamble = {
#    "font.family": "serif", # use serif/main font for text elements
#    "text.usetex": True,    # use inline math for ticks
#    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#    "pgf.texsystem": "xelatex",   # don't setup fonts from rc parameters
#    "pgf.preamble": " ".join([
#        r"\usepackage{units}",         # load additional packages
#        r"\usepackage{mathspec}",         # load additional packages
#        r"\setmainfont[" +\
#            "UprightFont = * ," +\
#            "ItalicFont = *Oblique," +\
#            "BoldFont = *Bold," +\
#            "Extension = .otf]{FreeSans}",
#        r"\setmathsfont(Digits,Latin,Greek)[" +\
#            "UprightFont = * ," +\
#            "ItalicFont = *Oblique," +\
#            "BoldFont = *Bold," +\
#            "Extension = .otf]{FreeSans}",
#        r"\setmathrm[" +\
#            "UprightFont = *," +\
#            "ItalicFont = *Oblique," +\
#            "BoldFont = *Bold," +\
#            "Extension = .otf]{FreeSans}",
#         ])
#}

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


the_ax = sns.relplot(
        data=subset_choosy_long.loc[subset_choosy_long["time"] < 20000]
        ,x="time"
        ,y="value"
        ,hue="variable"
        ,kind="line"
        )

the_ax.set(xlabel=r"Time, $t$",ylabel=r"Number infected, $I_{p\cdot}$ and $I_{c\cdot}$")
the_ax.savefig(fname="plot_choosy_vs_promiscuous.pdf", backend="pgf")
