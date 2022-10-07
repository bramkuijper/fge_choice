library("tidyverse")
library("khroma")

# overview plot of the various parameters
# - demographic feedbacks
# - cost of resistance
# - relative fecundity of good (but not best) plasmid

dat1 <- read_delim(file="summary_fge.csv",delim=";")
dat2 <- read_delim(file="summary_fge2.csv",delim=";")

dat <- bind_rows(dat1,dat2)

dat <- mutate(dat,
              Promiscuous=Ipg1+Ipg2
              ,Choosy=Icg1+Icg2
              )

dat <- filter(dat, demog_feedback == 0)

dat.pi.05 <- filter(dat, pi == 0.5)
dat.pi.1 <- filter(dat, pi == 1.0)
dat.pi.0.25 <- filter(dat, pi == 0.25)
dat.pi.0.05 <- filter(dat, pi == 0.05)

print(nrow(dat.pi.0.25))
print(nrow(dat.pi.0.05))



pg1 <- ggplot() +
  geom_contour(data=dat.pi.05
    ,mapping=aes(x=FG1, y=c, z=Choosy-Promiscuous)
               ,breaks=c(250), col="#000000") +
  geom_contour(data=dat.pi.1
    ,mapping=aes(x=FG1, y=c, z=Choosy-Promiscuous)
               ,breaks=c(250),col="#E69F00") +
  geom_contour(data=dat.pi.0.25
    ,mapping=aes(x=FG1, y=c, z=Choosy-Promiscuous)
               ,breaks=c(250),col="#56B4E9") +
  geom_contour(data=dat.pi.0.05
    ,mapping=aes(x=FG1, y=c, z=Choosy-Promiscuous)
               ,breaks=c(250),col="#009E73") +
  theme_classic() 
#  scale_fill_distiller(palette = "RdBu")
#  geom_raster(mapping=aes(fill=Choosy-Promiscuous)) +

ggsave(filename="overview_fge_choice.pdf")

dat.super <- read_delim(file="summary_superinfection.csv", delim=";")

dat.super <- mutate(dat.super
        ,Promiscuous=Sp + Ipg1 + Ipg2 + Ipg1g2
        ,Choosy=Sc + Icg1 + Icg2 + Icg1g2)

plot.super <- ggplot(data=dat.super
        ,mapping=aes(x=FG2/FG1,y=c)) +
        geom_tile(mapping=aes(fill=Choosy-Promiscuous)) +
        facet_grid(cols=vars(sigma))

ggsave(filename="overview_choice_superinfection.pdf")
