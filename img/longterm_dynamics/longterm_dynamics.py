#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import re, string
import os.path

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

panelcounter = 0

def single_panel(
        gsobject
        ,row
        ,col
        ,filename
        ,tmax
        ,ylim
        ,ylab=False
        ,yticklab=False
        ,xlab=False
        ,panel_label_loc="lower right"
        ,legend=False):

    global panelcounter
    
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

    colors = ["#831d50","#85cbee"]

    sns.set_palette(sns.color_palette(colors))
    
    plotje = sns.lineplot(
            data=subset_long.loc[subset_long["time"] < tmax]
            ,x="time"
            ,y="value"
            ,hue="variable"
            ,ax=the_ax
            ,legend=legend
            )

    sns.despine(ax=the_ax, top=True, right=True)

    the_ax.tick_params(bottom=True, 
            labelleft=yticklab,
            labelbottom=xlab)

    ylabel = ""

    if ylab:
        ylabel = r"Number infected, $I_{p\cdot}$ and $I_{c\cdot}$"

    xlabel = ""

    if xlab:
        xlabel=r"Time, $t$"


    the_ax.set(xlabel=xlabel
            ,ylabel=ylabel
            ,ylim=ylim)

    panel_label = AnchoredText(
            s=string.ascii_uppercase[panelcounter]
            ,frameon=False
            ,loc=panel_label_loc
)

    panelcounter += 1

    the_ax.add_artist(panel_label)

#    plt.close(plotje.fig)

    if legend:
        the_ax.legend(loc="lower center"
                ,fontsize=8
                ,frameon=False)

#    the_ax.savefig(fname="plot_choosy_vs_promiscuous.pdf")

def single_panel_types(
        gsobject
        ,row
        ,col
        ,filename
        ,tmax
        ,ylim
        ,ylab=False
        ,yticklab=False
        ,xlab=False
        ,panel_label_loc="lower right"
        ,legend=False):
    
    global panelcounter

    the_ax = plt.subplot(gsobject[row,col])

    # get line at which data file ends
    line_data_ends = find_parameter_line(
            filename = filename)

    # the data in which choosy outperforms anything else
    data = pd.read_csv(
            filepath_or_buffer=filename
            ,sep=";"
            ,nrows=line_data_ends - 1)

    data["G"] = data["Ipg1"] + data["Icg1"]
    data["B"] = data["Ipg2"] + data["Icg2"]

    # select columns we would like to plot
    subset = data[["time","G","B"]]

    # data currently in wide format, transform to long
    subset_long = pd.melt(frame=subset
            ,id_vars=["time"]
            ,value_vars=["G","B"])

    sns.set_palette(sns.color_palette("Set2"))
    
    plotje = sns.lineplot(
            data=subset_long.loc[subset_long["time"] < tmax]
            ,x="time"
            ,y="value"
            ,hue="variable"
            ,ax=the_ax
            ,legend=legend
            )

    sns.despine(ax=the_ax, top=True, right=True)

    the_ax.tick_params(bottom=True, 
            labelleft=yticklab,
            labelbottom=xlab)

    ylabel = ""

    if ylab:
        ylabel = r"Number infected with phage type $G$ vs $B$"

    xlabel = ""

    if xlab:
        xlabel=r"Time, $t$"


    the_ax.set(xlabel=xlabel
            ,ylabel=ylabel
            ,ylim=ylim)

    panel_label = AnchoredText(
            s=string.ascii_uppercase[panelcounter]
            ,frameon=False
            ,loc=panel_label_loc)

    panelcounter += 1

    the_ax.add_artist(panel_label)

    if legend:
        the_ax.legend(loc="lower center"
                ,fontsize=8
                ,frameon=False)

# loop through the different values of demog feedback
demog_feedback = [0,1]

file_dir = ["../growth_cp/no_demog_feedback/","../growth_cp"]

for demog_feedback_i in demog_feedback:

    # reset panel counter
    panelcounter = 0

    # make fig
    width = 10
    height = 10
    fig = plt.figure(figsize=(width,height))

    file_dir_i = file_dir[demog_feedback_i]

    file_no_infection = os.path.join(file_dir_i,"output_1")
    file_phageG = os.path.join(file_dir_i,"output_2")
    file_choosy_underperforms = os.path.join(file_dir_i,"output_3")
    file_choosy_outperforms = os.path.join(file_dir_i,"output_4")

    # start gridspec object
    gs = GridSpec(nrows=3, ncols=4,width_ratios=[1,1,0.25,1])

    ylim=[-30,1000]
    
    if demog_feedback_i == 0:
        ylim=[-30,1000]

    row_ctr = 0

    # first row: phage G only

    tmax=80000
    tmaxmax=1000000

    single_panel(gsobject=gs, 
            row=row_ctr,
            col=0,
            tmax=tmax,
            xlab=False,
            yticklab=True,
            ylim=ylim,
            filename=file_phageG)

    single_panel(gsobject=gs, 
            row=row_ctr,
            col=1,
            tmax=tmaxmax,
            xlab=False,
            ylim=ylim,
            panel_label_loc="upper right",
            filename=file_phageG)

    single_panel_types(gsobject=gs, 
            row=row_ctr,
            col=3,
            tmax=tmaxmax,
            xlab=False,
            yticklab=True,
            ylim=ylim,
            panel_label_loc="upper right",
            filename=file_phageG)

    row_ctr += 1

    # 2nd row: phage B only

    single_panel(gsobject=gs, 
            row=row_ctr,
            col=0,
            tmax=tmax,
            xlab=False,
            ylab=True,
            yticklab=True,
            legend=True,
            ylim=ylim,
            filename=file_choosy_underperforms)

    single_panel(gsobject=gs, 
            row=row_ctr,
            col=1,
            tmax=tmaxmax,
            xlab=False,
            ylim=ylim,
            filename=file_choosy_underperforms)

    single_panel_types(gsobject=gs, 
            row=row_ctr,
            col=3,
            tmax=tmaxmax,
            xlab=False,
            ylab=True,
            yticklab=True,
            legend=True,
            ylim=ylim,
            filename=file_choosy_underperforms)

    row_ctr += 1

    # row 3: mixture

    single_panel(gsobject=gs, 
            row=row_ctr,
            col=0,
            tmax=tmax,
            xlab=True,
            yticklab=True,
            ylim=ylim,
            filename=file_choosy_outperforms)

    single_panel(gsobject=gs, 
            row=row_ctr,
            col=1,
            tmax=tmaxmax,
            xlab=True,
            ylim=ylim,
            filename=file_choosy_outperforms)

    single_panel_types(gsobject=gs, 
            row=row_ctr,
            col=3,
            tmax=tmaxmax,
            xlab=True,
            yticklab=True,
            ylim=ylim,
            filename=file_choosy_outperforms)


    plt.figtext(x=-0.05, y= 0.8, s= "Only MGE B\n" + "→")
    plt.figtext(x=-0.05, y= 0.5, s= "Only MGE G\n" + "→")
    plt.figtext(x=-0.05, y= 0.15, s= "Both \n" + "MGEs" + "→")

    figname = "plot_longterm_choosy_vs_promiscuous_dg" + str(demog_feedback_i) + ".pdf"

    fig.savefig(
            fname=figname
            ,bbox_inches="tight"
            )
