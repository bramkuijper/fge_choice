#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import re, string

from matplotlib import rcParams
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import AnchoredText
import matplotlib as mpl
import matplotlib.pyplot as plt

fontpath = "/System/Library/Fonts/Supplemental/" 
#fontpath = "/usr/local/texlive/2022/texmf-dist/fonts/opentype/public/stix2-otf/"

# set basic theme
print(mpl.matplotlib_fname())
mpl.rcParams.update(mpl.rcParamsDefault)


sns.set_theme(style="ticks")
mpl.rcParams["mathtext.fontset"] = "dejavusans"
#
## see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
#pgf_with_custom_preamble = {
#    "font.family": "sans-serif", # use serif/main font for text elements
#    "pgf.texsystem": "xelatex",
#    "text.usetex": True,    # use inline math for ticks
#    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#    "pgf.preamble": "\n".join([
#        r"\usepackage{units}",         # load additional packages
#        r"\usepackage{mathspec}",         # load additional packages
#        r"\setmainfont[" +\
#            "Path = " + fontpath + "," +\
#            "UprightFont = *-Regular," +\
#            "ItalicFont = *-Italic," +\
#            "BoldFont = *-Bold," +\
#            "Extension = .otf]{STIXTwoText}",
#        r"\setmathfont[" +\
#            "Path = " + fontpath + "," +\
#            "UprightFont = *-Regular," +\
#            "ItalicFont = *-Italic," +\
#            "BoldFont = *-Bold," +\
#            "Extension = .otf]{STIXTwoText}",
#        r"\setmathrm[" +\
#            "Path = " + fontpath + "," +\
#            "UprightFont = *-Regular," +\
#            "ItalicFont = *-Italic," +\
#            "BoldFont = *-Bold," +\
#            "Extension = .otf]{STIXTwoText}",
#         ])
#}
#
#mpl.rcParams.update(pgf_with_custom_preamble)


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



def single_panel(
        gsobject
        ,row
        ,col
        ,filename
        ,legend=False):
    
    the_ax = plt.subplot(gsobject[row,col])

    # get line at which data file ends
    line_data_ends = find_parameter_line(
            filename = filename)

    # the data in which choosy outperforms anything else
    data = pd.read_csv(
            filepath_or_buffer=filename
            ,sep=";"
            ,nrows=line_data_ends - 1)

    data["Promiscuous"] = data["Ipg1"] + data["Ipg2"]
    data["Choosy"] = data["Icg2"] + data["Icg1"]

    # select columns we would like to plot
    subset = data[["time","Choosy","Promiscuous"]]

    # data currently in wide format, transform to long
    subset_long = pd.melt(frame=subset
            ,id_vars=["time"]
            ,value_vars=["Choosy","Promiscuous"])

    tmax = 1800

    colors = ["#cc79a7","#3873b2"]

    sns.set_palette(sns.color_palette(colors))
    

    plotje =sns.lineplot(
            data=subset_long.loc[subset_long["time"] < tmax]
            ,x="time"
            ,y="value"
            ,hue="variable"
            ,ax=the_ax
            ,legend=legend
            )

    remove_y = col != 0

    sns.despine(ax=the_ax, left=remove_y, top=True, right=True)

    the_ax.tick_params(left=not remove_y, 
            labelleft=not remove_y)

    ylab = ""

    if not remove_y:
        ylab = r"Number infected, $I_{p\cdot}$ and $I_{c\cdot}$"

    the_ax.set(xlabel=r"Time, $t$"
            ,ylabel=ylab)

    panel_label = AnchoredText(
            s=string.ascii_uppercase[col]
            ,frameon=False
            ,loc="lower right")

    the_ax.add_artist(panel_label)

#    plt.close(plotje.fig)

    if legend:
        the_ax.legend(loc="center right"
                ,fontsize=8
                ,frameon=False)

#    the_ax.savefig(fname="plot_choosy_vs_promiscuous.pdf")


# file name 
file_choosy_outperforms = "choosy_outperforms_short_term.csv"
file_choosy_underperforms = "choosy_underperforms.csv"

width = 10
height = 3
fig = plt.figure(figsize=(width,height))

# start gridspec object
gs = GridSpec(nrows=1, ncols=4)

single_panel(gsobject=gs, 
        row=0,
        col=0,
        filename=file_choosy_underperforms)

single_panel(gsobject=gs, 
        row=0,
        col=1,
        filename=file_choosy_outperforms,
        legend=True)

fig.savefig(
        fname="plot_choosy_vs_promiscuous.pdf"
        ,bbox_inches="tight"
        )
