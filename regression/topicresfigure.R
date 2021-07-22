#topicmodel_figure


topic.res.df <- read_csv('./regression/topicmodel_deltas5.csv')

topic.cats <- read.csv('./topicdata/k200standardcoding.tsv', sep = '\t')

library(ggbeeswarm)

topic.res.df %>%
  mutate(topic.label = topic.cats$label,
         category = topic.cats$category) %>%
  # filter(category %in% 
  #          c("dialect / language", "event", "genre", 
  #            "human institutions, practices, or relationships",
  #            "nationalities, regions, or ethnicities",
  #            "physical description", "technology",
  #            "uncategorized dimension of style")) %>%
  ggplot(aes(x = reorder(category, adjdelta), y = adjdelta, fill = category,
             label = topic.label)) + 
  geom_quasirandom(shape = 21) + 
  geom_boxplot(alpha = .5) + 
  # geom_text_repel(seed = 12345, size = 2, max.iter = 10000) +
  # facet_grid(category~., scales = "free_y", space = "free",
  #            switch = "y") + 
  coord_flip() + 
  labs(x = "LIWC Category", y = "Adjusted Delta") + 
  theme_bw() + 
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle = 0)) + 
  labs(x = "")







