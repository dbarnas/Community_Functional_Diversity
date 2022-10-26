### Get color palette from Mo'orea using
### LTER Palette Finder
### https://github.com/lter/lterpalettefinder#readme

#devtools::install_github("lter/lterpalettefinder")

library(lterpalettefinder)
library(here)

palette_find()

sunsetPal <- palette_extract(image = here("Data","Images","MooreaSunset1.jpg"))

sortedPal <- palette_sort(sunsetPal)

palette_demo(sortedPal)

palette_ggdemo(sortedPal)
