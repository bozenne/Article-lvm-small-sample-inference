source("ANALYSIS.R")

dtL.vitamin[, group := factor(grp, levels = c("C","T"), labels = c("Control","Treatment"))]
gg.spaguetti <- ggplot(dtL.vitamin, aes(x = as.factor(week), y = weight0, group = animal,
                                        color = animal))
gg.spaguetti <- gg.spaguetti + geom_line(size = 2) + geom_point(size = 3)
gg.spaguetti <- gg.spaguetti + facet_grid(~group, labeller = label_both)
gg.spaguetti <- gg.spaguetti + xlab("week") + ylab("weight")
gg.spaguetti <- gg.spaguetti + theme(text = element_text(size=20))

ggsave(gg.spaguetti, filename = file.path("Figures","fig-spaguetti.pdf"))
