library("tidyverse")
library("khroma")

# overview plot of the various parameters
# - demographic feedbacks
# - cost of resistance
# - relative fecundity of good (but not best) plasmid

dat <- read_delim(file="summary_fge.csv",delim=";")

dat <- mutate(dat,
              Promiscuous=Ipg1+Ipg2
              ,Choosy=Icg1+Icg2
              )

dat <- filter(dat, demog_feedback == 0)

dat.pi.05 <- filter(dat, pi == 0.5)
dat.pi.1 <- filter(dat, pi == 1.0)



pg1 <- ggplot() +
  geom_contour(data=dat.pi.05
    ,mapping=aes(x=FG1, y=c, z=Choosy-Promiscuous)
               ,breaks=c(250)) +
  geom_contour(data=dat.pi.1
    ,mapping=aes(x=FG1, y=c, z=Choosy-Promiscuous)
               ,breaks=c(250),col="red") +
  theme_classic() 
#  scale_fill_distiller(palette = "RdBu")
#  geom_raster(mapping=aes(fill=Choosy-Promiscuous)) +

ggsave(filename="overview_fge_choice.pdf")
