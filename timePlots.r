library(dplyr)
library(ggplot2)
library(reshape2)

dat <- read.csv("~/lauraMouse/timeplots/allTPdata.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control",expression(italic("Villin-cre, Bai1"^{"WT/Tg"})))
#myLabels <- c("Control","Bai1")
myLevels <- c("WT","Tg")
myColors <- c("yellow","cyan")
OUTDIR <- "~/lauraMouse/timeplots/"

dat <- read.csv("~/lauraMouse/samantha/villinCreElmo20180716.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control",expression(paste(italic("Villin-cre, ELMO1"^{"FL/FL"}),italic("ELMO2"^{"FL/FL"}),sep=" ")))
#myLabels <- c("Control","Bai1")
myLevels <- c("WT","KO")
myColors <- c("yellow","purple")
OUTDIR <- "~/lauraMouse/samantha/"

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
		   eduDiffSE = sqrt(eduSE.uv^2 + eduSE.nouv^2),
		   caspDiff = casp.uv - casp.nouv,
		   caspDiffSE = sqrt(caspSE.uv^2 + caspSE.nouv^2),
		   tunelDiff = tunel.uv - tunel.nouv,
		   tunelDiffSE = sqrt(tunelSE.uv^2 + tunelSE.nouv^2),
		   eduAreaDiff = eduArea.uv - eduArea.nouv, 
		   eduAreaDiffSE = sqrt(eduAreaSE.uv^2 + eduAreaSE.nouv^2),
		   caspAreaDiff = caspArea.uv - caspArea.nouv,
		   caspAreaDiffSE = sqrt(caspAreaSE.uv^2 + caspAreaSE.nouv^2),
		   tunelAreaDiff = tunelArea.uv - tunelArea.nouv,
		   tunelAreaDiffSE = sqrt(tunelAreaSE.uv^2 + tunelAreaSE.nouv^2)
		   ) %>%
	select(mouse,genotype,TIMEPOINT,eduDiff,caspDiff,tunelDiff,eduAreaDiff,caspAreaDiff,tunelAreaDiff,
		   eduDiffSE, caspDiffSE, tunelDiffSE, eduAreaDiffSE, caspAreaDiffSE, tunelAreaDiffSE)

diffGraph <- melt(diffs,id.vars=c("mouse","genotype","TIMEPOINT"))
diffGraph$variable <- factor(diffGraph$variable, levels=c("eduDiff","caspDiff","tunelDiff","eduAreaDiff","caspAreaDiff","tunelAreaDiff"),
                                 labels=c("EdU","Caspase-3","TUNEL","EdU Area","Caspase-3 Area","TUNEL Area"))


doTimePlots(diffGraph, "EdU",       "UV Intensity - NOUV Intensity", "EDU",       "eduIntensity.png")
doTimePlots(diffGraph, "Caspase-3", "UV Intensity - NOUV Intensity", "Caspase-3", "caspaseIntensity.png")
doTimePlots(diffGraph, "TUNEL",     "UV Intensity - NOUV Intensity", "TUNEL",     "tunelIntensity.png")

doTimePlots(diffGraph, "EdU Area",       "UV Area - NOUV Area", "EDU",       "eduArea.png")
doTimePlots(diffGraph, "Caspase-3 Area", "UV Area - NOUV Area", "Caspase-3", "caspaseArea.png")
doTimePlots(diffGraph, "TUNEL Area",     "UV Area - NOUV Area", "TUNEL",     "tunelArea.png")

# Plots for each mouse
controlDiffs <- filter(diffs, genotype == "Control")
genoDiffs <- filter(diffs, genotype != "Control")

doMousePlots(controlDiffs, "UV Intensity - NOUV Intensity", "Caspase-3 - Control", "caspDiff", "caspDiffSE", "caspaseIntensityAllMice-Control.png")
doMousePlots(controlDiffs, "UV Intensity - NOUV Intensity", "EdU - Control", "eduDiff", "eduDiffSE", "eduIntensityAllMice-Control.png")
doMousePlots(controlDiffs, "UV Intensity - NOUV Intensity", "TUNEL - Control", "tunelDiff", "tunelDiffSE", "tunelIntensityAllMice-Control.png")

doMousePlots(genoDiffs, "UV Intensity - NOUV Intensity", "Caspase-3 - Bai1", "caspDiff", "caspDiffSE", "caspaseIntensityAllMice-Bai1.png")
doMousePlots(genoDiffs, "UV Intensity - NOUV Intensity", "EdU - Bai1", "eduDiff", "eduDiffSE", "eduIntensityAllMice-Bai1.png")
doMousePlots(genoDiffs, "UV Intensity - NOUV Intensity", "TUNEL - Bai1", "tunelDiff", "tunelDiffSE", "tunelIntensityAllMice-Bai1.png")

doMousePlots(controlDiffs, "UV Area - NOUV Area", "Caspase-3 - Control", "caspAreaDiff", "caspAreaDiffSE", "caspaseAreaAllMice-Control.png")
doMousePlots(controlDiffs, "UV Area - NOUV Area", "EdU - Control", "eduAreaDiff", "eduAreaDiffSE", "eduAreaAllMice-Control.png")
doMousePlots(controlDiffs, "UV Area - NOUV Area", "TUNEL - Control", "tunelAreaDiff", "tunelAreaDiffSE", "tunelAreaAllMice-Control.png")

doMousePlots(genoDiffs, "UV Area - NOUV Area", "Caspase-3 - Bai1", "caspAreaDiff", "caspAreaDiffSE", "caspaseAreaAllMice-Bai1.png")
doMousePlots(genoDiffs, "UV Area - NOUV Area", "EdU - Bai1", "eduAreaDiff", "eduAreaDiffSE", "eduAreaAllMice-Bai1.png")
doMousePlots(genoDiffs, "UV Area - NOUV Area", "TUNEL - Bai1", "tunelAreaDiff", "tunelAreaDiffSE", "tunelAreaAllMice-Bai1.png")


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

doMousePlots <- function(myData, myYlabel, myTitle, colName, errName, fileName) {
	myDodge <- 0.3
    ggplot() +
	  geom_hline(aes(yintercept=0), linetype="dotted") +
	  geom_line(aes(x = TIMEPOINT, y = get(colName), group=mouse, colour=mouse), size = 1.5, linetype = 2, position=position_dodge(myDodge), data=myData) +
      geom_errorbar(aes(x = TIMEPOINT, fill=mouse, colour=mouse, ymin=get(colName)-get(errName), ymax=get(colName)), width=.2, position=position_dodge(myDodge), data=myData) +
      geom_errorbar(aes(x = TIMEPOINT, fill=mouse, colour=mouse, ymin=get(colName), ymax=get(colName)+get(errName)), width=.2, position=position_dodge(myDodge), data=myData) +
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
	  scale_x_discrete(limits=c("3HR","6HR","12HR","24HR"))
      ggsave(filename=fileName, path = OUTDIR, width=10, height=7, bg="transparent")
}
