library(dplyr)
library(ggplot2)
library(reshape2)

dat <- read.csv("~/lauraMouse/timeplots/allTPdata.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control",expression(italic("Villin-cre, Bai1"^{"WT/Tg"})))
#myLabels <- c("Control","Bai1")
myLevels <- c("WT","Tg")
myColors <- c("yellow","cyan")
OUTDIR <- "~/lauraMouse/timeplots2/"

myData <- dat %>%
  mutate(genotype = factor(GENOTYPE, levels = myLevels, labels=myLabels)) %>%
  rename(mouse = MOUSE.., treatment = TREATMENT) %>%
#  filter(mouse != "BM30") %>%
  group_by(mouse,treatment,genotype,TIMEPOINT) %>%
  summarise( count = n(),
		   area = mean(AreaShape_Area),
		   areaSE = sd(AreaShape_Area)/sqrt(count),
		   perimeter = mean(AreaShape_Perimeter),
		   perimeterSE = sd(AreaShape_Perimeter)/sqrt(count),
		   casp = mean(CASP.intensity.corrected),
		   caspSE = sd(CASP.intensity.corrected)/sqrt(count),
		   edu = mean(EDU.Intensity.corrected),
		   eduSE = sd(EDU.Intensity.corrected)/sqrt(count),
		   tunel = mean(TUNEL.Intensity.Corrected),
		   tunelSE = sd(TUNEL.Intensity.Corrected)/sqrt(count),
		   caspArea = mean(CASP.area.corrected),
		   caspAreaSE = sd(CASP.area.corrected)/sqrt(count),
		   eduArea = mean(EDU.area.corrected),
		   eduAreaSE = sd(EDU.area.corrected)/sqrt(count),
		   tunelArea = mean(TUNEL.area.Corrected),
		   tunelAreaSE = sd(TUNEL.area.Corrected)/sqrt(count) )

uv <- filter(myData, treatment == "UV")
nouv <- filter(myData, treatment == "NOUV")
diffs <- inner_join(uv, nouv, by=c("mouse","genotype","TIMEPOINT"), suffix=c(".uv",".nouv")) %>%
	ungroup() %>%
	mutate(eduDiff = edu.uv - edu.nouv,
		   caspDiff = casp.uv - casp.nouv,
		   tunelDiff = tunel.uv - tunel.nouv,
		   eduAreaDiff = eduArea.uv - eduArea.nouv, 
		   caspAreaDiff = caspArea.uv - caspArea.nouv,
		   tunelAreaDiff = tunelArea.uv - tunelArea.nouv
		   ) %>%
	select(mouse,genotype,TIMEPOINT,eduDiff,caspDiff,tunelDiff,eduAreaDiff,caspAreaDiff,tunelAreaDiff)

diffGraph <- melt(diffs,id.vars=c("mouse","genotype","TIMEPOINT"))
diffGraph$variable <- factor(diffGraph$variable, levels=c("eduDiff","caspDiff","tunelDiff","eduAreaDiff","caspAreaDiff","tunelAreaDiff"),
                                 labels=c("EdU","Caspase-3","TUNEL","EdU Area","Caspase-3 Area","TUNEL Area"))


doTimePlots(diffGraph, "EdU",       "UV Intensity - NOUV Intensity", "EDU",       "eduIntensity.png")
doTimePlots(diffGraph, "Caspase-3", "UV Intensity - NOUV Intensity", "Caspase-3", "caspaseIntensity.png")
doTimePlots(diffGraph, "TUNEL",     "UV Intensity - NOUV Intensity", "TUNEL",     "tunelIntensity.png")

doTimePlots(diffGraph, "EdU Area",       "UV Area - NOUV Area", "EDU",       "eduArea.png")
doTimePlots(diffGraph, "Caspase-3 Area", "UV Area - NOUV Area", "Caspase-3", "caspaseArea.png")
doTimePlots(diffGraph, "TUNEL Area",     "UV Area - NOUV Area", "TUNEL",     "tunelArea.png")

# Plot all the data
for (i in c("casp","caspArea","edu","eduArea","tunel","tunelArea")) {
	for (t in c("3HR","6HR","12HR","24HR")) {
		comparePlotsByTime(myData, i, t)
	}
}

# Do the anova's and write out the tukey's
# TODO: Change the label to make the csv readable
pVals <- data.frame("label"=c("genotype","TIMEPOINT","interaction"))
writeHeader <- TRUE
for (x in c("eduDiff","caspDiff","tunelDiff","eduAreaDiff","caspAreaDiff","tunelAreaDiff")) {
  a <- aov(get(x) ~ genotype*TIMEPOINT, data=diffs)
  pVals[x] <- summary(a)[[1]][["Pr(>F)"]][1:3]
  tky = cbind(x, TukeyHSD(a)$`genotype:TIMEPOINT`)
  write.table(tky, file=paste(OUTDIR,"/tukey.csv", sep=""), append=!writeHeader, col.names=writeHeader, quote=FALSE, sep=",")
  writeHeader <- FALSE
}


doTimePlots <- function(allDiffData,myVar,myYlabel,myTitle,fileName,myDodge=0.4) {

	diffGraph <- filter(allDiffData, variable == myVar)

    diffError <- diffGraph %>%
      group_by(genotype,TIMEPOINT) %>%
      summarise( m = mean(value), se = sd(value)/sqrt(n()) )

    ggplot() + 
	  geom_hline(aes(yintercept=0), linetype="dotted") +
	  geom_line(aes(x = TIMEPOINT, y = m, group=genotype, colour=genotype), size = 1.5, linetype = 2, position=position_dodge(myDodge), data=diffError) +
      geom_dotplot(aes(y = value, x = TIMEPOINT, fill = genotype), data = diffGraph,
                   binaxis='y', stackdir='center', dotsize=1, stackgroups = TRUE, position=position_dodge(myDodge)) +
      geom_errorbar(aes(x = TIMEPOINT, fill=genotype, ymin=m-se, ymax=m), width=.2, position=position_dodge(myDodge), data=diffError) +
      geom_errorbar(aes(x = TIMEPOINT, fill=genotype, ymin=m, ymax=m+se), width=.2, position=position_dodge(myDodge), data=diffError) +
      theme_classic(base_size=14) +
      theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            legend.text.align = 0,
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
            rect = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent", colour= NA),
            panel.background = element_rect(fill = "transparent", colour= NA),
            legend.text=element_text(size=14),
			plot.title = element_text(hjust = 0.5)
		   ) +
      labs(y=myYlabel,x="") +
	  ggtitle(myTitle) +
	  scale_x_discrete(limits=c("3HR","6HR","12HR","24HR")) +
      scale_fill_manual(values=myColors, name="", labels=myLabels) +
      scale_colour_manual(values=myColors, name="", labels=myLabels)
      ggsave(filename=fileName, path = OUTDIR, width=10, height=7, bg="transparent")
}


#allData <- myData
#colName <- "casp"
#tp <- "3HR"

comparePlotsByTime <- function(allData, colName, tp) {
 compData <- myData %>% filter(TIMEPOINT == tp)
 fileName <- paste(colName,tp,"png",sep=".")
 ggplot(compData, aes(x = treatment, y = get(colName), label=mouse, colour = genotype, group = mouse, fill = genotype)) +
   geom_line(size=1.5) +
   geom_text(colour="black") +
   #geom_point(colour="black",pch=21, size=5) +
   facet_wrap(~ genotype, labeller=label_parsed) +
   scale_colour_manual(name="", values=myColors, guide=FALSE) +
   scale_fill_manual(name="", values=myColors, guide=FALSE) +
    theme_classic(base_size=14) +
    theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.text.align = 0,
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          rect = element_rect(fill = "transparent"),
		  strip.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour= NA),
          panel.background = element_rect(fill = "transparent", colour= NA),
          legend.text=element_text(size=14)) +
    labs(y=paste(colName,tp),x="Treatment")
  ggsave(paste(OUTDIR, fileName, sep="/"), bg="transparent")
}
