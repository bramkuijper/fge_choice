library("tidyverse")
library("khroma")

# overview plot of the various parameters
# - demographic feedbacks
# - strength of resistance
# - loss rate
# all in the context of superinfection, as higher loss rate may result in things transmitting through the system prolongedly, as superinfecteds cannot act as a reservoir


dat.super <- read_delim(file="summary_gamma_pi.csv", delim=";")

dat.super <- mutate(dat.super
        ,Promiscuous=Sp + Ipg1 + Ipg2 + Ipg1g2
        ,Choosy=Sc + Icg1 + Icg2 + Icg1g2)

plot.super <- ggplot(data=dat.super) +
        geom_contour(
        mapping=aes(x=gammaCG1,y=sigma, z=Choosy-Promiscuous)
                ,breaks=c(0)
                ) +
        facet_grid(cols=vars(demog_feedback), rows=vars(pi))

ggsave(filename="overview_choice_superinfection.pdf")


dat.super.2 <- read_delim(file="summary_gamma_cost.csv", delim=";")

dat.super.2 <- mutate(dat.super.2
        ,Promiscuous=Sp + Ipg1 + Ipg2 + Ipg1g2
        ,Choosy=Sc + Icg1 + Icg2 + Icg1g2
        ,Diff=ifelse(Choosy - Promiscuous > 0,1,-1))

smooth_subset <- function(the_data)
{
  # use loess smoothers
  #  https://stackoverflow.com/questions/9874973/smoothness-of-contour-lines-in-contourplot
  b <- loess(Diff ~ gammaCG1 + sigma, data = the_data,span=0.03)
  
  smooth <- predict(b, 
          newdata=data.frame(
                  gammaCG1=the_data$gammaCG1
                  ,sigma=the_data$sigma))
  
  the_data$smooth <- smooth 
  return(the_data)
}

all_sub_dat <- NULL

for (c.i in sort(unique(dat.super.2$c)))
{
  sub_dat <- filter(.data=dat.super.2,c == c.i)
  
  sub_dat <- smooth_subset(the_data=sub_dat)
  
  all_sub_dat <- bind_rows(all_sub_dat,sub_dat)
}

plot.super.2 <- ggplot(data=all_sub_dat) +
        geom_contour(data=filter(all_sub_dat,c==0.01)
                ,mapping=aes(x=gammaCG1,y=sigma,z=smooth)
                ,breaks=c(0)
                ,colour="#66c2a5") +
        geom_contour(data=filter(all_sub_dat,c==0.05)
                ,mapping=aes(x=gammaCG1,y=sigma,z=smooth)
                ,breaks=c(0)
                ,colour="#fc8d62") +
        geom_contour(data=filter(all_sub_dat,c==0.1)
                ,mapping=aes(x=gammaCG1,y=sigma,z=smooth)
                ,breaks=c(0)
                ,colour="#8da0cb") +
        geom_contour(data=filter(all_sub_dat,c==0.2)
                ,mapping=aes(x=gammaCG1,y=sigma,z=smooth)
                ,breaks=c(0)
                ,colour="#e78ac3") +
        theme_classic() +
        xlab("FGE loss rate, \u03b3") + 
        ylab("Probability of superinfection, \u03c3") +
  annotate("text", x=4.5, y=0.018, label="c = 0.5") + 
  annotate("text", x=3.05,y=0.085,label="c = 0.01", angle=60) + 
  annotate("text", x=3.6,y=0.08,label="c = 0.05",angle=50) +
  annotate("text", x=4,y=0.062,label="c = 0.1",angle=30) +
  annotate("text", x=4.5,y=0.005,label="Choice") +
  annotate("text", x=1.5,y=0.075,label="No choice") 

ggsave(filename="overview_choice_superinfection_vary_cost.pdf"
       ,device=cairo_pdf
       ,width=5
       ,height=5)
