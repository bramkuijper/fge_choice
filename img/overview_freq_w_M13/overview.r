library("tidyverse")

# read in the file
f <- read_delim(file="summary.csv",delim = ";")

f_m <- mutate(f
              ,Gdelta = pGc - pGp
              ,Bdelta = pBc - pBp
              ,time = ifelse(time > 1000,10000,time))

# take mixed infections (i.e., those where G1 and G2 have nonzero
# frequencies at t = 0)
mixed_only <- f_m %>% filter(IPG2 > 0 & IPG1 > 0 & time > 0 & time < 10000)


mixed_only_l <- pivot_longer(data=mixed_only
                     ,cols = c("Gdelta","Bdelta")
                     ,values_to = "conditional_value"
                     ,names_to = "conditional")

ggplot(data = mixed_only_l
       ,mapping = aes(x=FG1,y=dSP)) +
    geom_tile(mapping = aes(fill=conditional_value)) +
    facet_grid(cols=vars(conditional),rows=vars(time)) +
    theme_classic(base_size=16) +
    scale_fill_distiller(palette = "Spectral")

ggsave("overview_plot_conditionals.pdf")
    


