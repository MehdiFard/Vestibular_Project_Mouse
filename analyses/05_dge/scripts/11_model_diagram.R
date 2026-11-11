# 11 - A diagram of what each model formula actually compares.
#
# The formula ~ strain + batch + age lists PREDICTORS. It does not say what is
# being compared - that is decided separately, by the contrast you ask for.
# Everything else in the formula is there to be held constant.

source("../../analyses/05_dge/scripts/00_config.R")

library(data.table)
library(ggplot2)

# One strain's four libraries: two organs x two ages.
cells <- CJ(organ = c("Saccule", "Utricle"), age = c("w8", "w22"))
cells[, x := ifelse(organ == "Saccule", 1, 2)]
cells[, y := ifelse(age == "w8", 2, 1)]
cells[, lab := paste0(substr(organ, 1, 3), "\n", age)]

panel <- function(arrows, title, formula, question, compared) {
  d <- copy(cells)[, panel := title]
  a <- copy(arrows)[, panel := title]
  list(cells = d, arrows = a,
       meta = data.table(panel = title, formula = formula,
                         question = question, compared = compared))
}

# The organ contrast: compare left to right, inside one age band.
arr_organ <- data.table(x = c(1, 1), xend = c(2, 2), y = c(2, 1), yend = c(2, 1))
# The age contrast: compare bottom to top, inside one organ column.
arr_age   <- data.table(x = c(1, 2), xend = c(1, 2), y = c(2, 2), yend = c(1, 1))

# Pull each arrow back from the box edges, along its own direction.
shrink <- function(a, dx = 0.38, dy = 0.33) {
  a[, `:=`(x2    = x    + sign(xend - x) * dx,
           xend2 = xend - sign(xend - x) * dx,
           y2    = y    + sign(yend - y) * dy,
           yend2 = yend - sign(yend - y) * dy)]
  a[]
}
arr_organ <- shrink(arr_organ)
arr_age   <- shrink(arr_age)

p1 <- panel(arr_organ, "A. Organ contrast",
            "~ strain + organ",
            "Within one animal, does the utricle differ from its own saccule?",
            "organ")
p2 <- panel(arr_age, "B. Age contrast",
            "~ strain + batch + age",
            "Within one strain, does 22 weeks differ from 8 weeks?",
            "age")

all_cells  <- rbindlist(list(p1$cells, p2$cells))
all_arrows <- rbindlist(list(p1$arrows, p2$arrows))
meta       <- rbindlist(list(p1$meta, p2$meta))

all_cells[, panel := factor(panel, levels = meta$panel)]
all_arrows[, panel := factor(panel, levels = meta$panel)]

p <- ggplot(all_cells, aes(x, y)) +
  geom_tile(aes(fill = organ), width = 0.72, height = 0.62,
            colour = "grey30", alpha = 0.55) +
  geom_segment(data = all_arrows, aes(x = x2, xend = xend2, y = y2, yend = yend2),
               arrow = arrow(length = unit(0.22, "cm"), ends = "both"),
               linewidth = 0.9, colour = "#B2182B", inherit.aes = FALSE) +
  geom_text(aes(label = lab), size = 3.5, lineheight = 0.95) +
  facet_wrap(~ panel) +
  scale_fill_manual(values = ORGAN_COLS, guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("Saccule", "Utricle"), limits = c(0.5, 2.5)) +
  scale_y_continuous(breaks = 1:2, labels = c("22 weeks", "8 weeks"), limits = c(0.5, 2.5)) +
  labs(title = "The formula lists what is held constant. The contrast picks what is compared.",
       subtitle = paste("Each box is one library from a single strain. Red arrows are the comparison",
                        "the contrast asks for."),
       x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11))

ggsave(file.path(FIGS, "29_model_diagram.png"), p,
       width = 10, height = 4.6, dpi = 300, bg = "white")

# ---------------------------------------------------------------- companion table
tab <- data.table(
  Model = c("~ strain + organ", "~ batch + age", "~ strain + batch + age",
            "~ batch + organ*age"),
  `Contrast asked for` = c("organ: Utricle vs Saccule", "age: w22 vs w8",
                           "age: w22 vs w8", "organ:age interaction"),
  `Held constant` = c("strain (and batch, which travels with it)",
                      "batch", "strain and batch", "batch"),
  `Question` = c("Does the utricle differ from the saccule, within an animal?",
                 "Does 22 weeks differ from 8 weeks, across animals?",
                 "Does 22 weeks differ from 8 weeks, within a strain?",
                 "Does the organ difference itself change with age?"))
fwrite(tab, file.path(TABLES, "model_formula_guide.tsv"), sep = "\t")
print(tab)

cat("\nDone.\n")
