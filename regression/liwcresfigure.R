library(ggrepel)


#Load in LIWC and cross-validated results
res.df <- read_csv('crossvalidated_delta_knit.tsv')
d <- read_tsv("../liwc/liwc_w_meta.tsv")

#Making category labels for presentation
liwc_categories <- tibble(liwc_var = names(d)[9:89]) %>%
  mutate(category = c("word count", rep("summary language", 7),
                      rep("linguistic dimensions", 15),
                      rep("other grammar", 6), 
                      rep("affective processes", 6),
                      rep("social processes", 5), 
                      rep("cognitive processes", 7),
                      rep("perceptual processes", 4),
                      rep("biological processes", 5),
                      rep("drives", 6),
                      rep("time orientations", 4),
                      rep("relativity", 4),
                      rep("personal concerns", 6),
                      rep("informal language", 5)))
rm(d) # Don't need this anymore

#

pdf("./regression/liwcresplot.pdf") 
# 2. Create a plot

res.df %>%
  mutate(xmin = adjdelta - .1, xmax = adjdelta + .1) %>%
  left_join(liwc_categories, by = c("depvar"="liwc_var")) %>%
  ggplot(aes(x = reorder(category, adjdelta), y = adjdelta, fill = category,
             label = depvar)) + 
  geom_point(shape = 21) + 
  geom_text_repel(seed = 12345, size = 3, max.iter = 10000,
                  nudge_x = .3) + 
  facet_grid(category~., scales = "free_y", switch = "y") + 
  coord_flip() + 
  labs(x = "LIWC Category", y = "Adjusted Delta") + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle = 0),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

dev.off() 

