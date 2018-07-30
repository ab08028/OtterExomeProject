##### Faststructure outputs SVG plots 
# download them and copy to : /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/FASTSTRUCTURE
wd="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/plots/FASTSTRUCTURE"
# install.packages("rsvg")
# install.packages("svglite")
require(rsvg)
require(svglite)
test <- "/Users/annabelbeichman/Downloads/XX.test.faststructure_plot.svg"
test1 <- rsvg(test)
png::writePNG(test1, "test1.png")
