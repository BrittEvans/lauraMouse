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

dat <- read.csv("~/lauraMouse/data/irradiation.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplots/"

dat <- read.csv("~/lauraMouse/data/CsIrradiatedNov2018.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplotsLaura/"

dat <- read.csv("~/lauraMouse/data/irradiationCombined.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplotsCombined/"

dat <- read.csv("~/lauraMouse/data/CsIrradiatedJan2019.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplotsLauraJan2019/"

dat <- read.csv("~/lauraMouse/data/CsIrradiatedNovAndJan2019.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplotsLauraNovAndJan2019/"

dat <- read.csv("~/lauraMouse/data/cx18cx22.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplotsCX18CX22/"

dat <- read.csv("~/lauraMouse/data/CsIrradiatedApril2019.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control","700 rads")
myLevels <- c("NOIR","IR")
myColors <- c("yellow","black")
OUTDIR <- "~/lauraMouse/timeplotsApril2019/"


dat2 <- read.csv("~/lauraMouse/data/CsIrradiatedNovAndJan2019.csv")
OUTDIR <- "~/lauraMouse/timeplotsAllCombined/"

dat1 <- dat %>%
  filter(Metadata_Time != "3HR") %>%
  filter(Metadata_Well %in% c("1","2"))

dat <- bind_rows(dat1, dat2)
names(dat1)
names(dat2)



myData <- dat %>%
  mutate(treatment = factor(Metadata_Treatment, levels=myLevels, labels=myLabels)) %>%
  rename(mouse = Metadata_Mouse, genotype = GENOTYPE) %>%
  group_by(mouse,treatment,Metadata_Time) %>%
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

uv <- filter(myData, treatment == "700 rads")
nouv <- filter(myData, treatment == "Control")
diffs <- inner_join(uv, nouv, by=c("mouse","Metadata_Time"), suffix=c(".uv",".nouv")) %>%
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
		   tunelAreaDiffSE = sqrt(tunelAreaSE.uv^2 + tunelAreaSE.nouv^2),
		   areaDiff = area.uv - area.nouv,
		   areaDiffSE = sqrt(areaSE.uv^2 + areaSE.nouv^2),
		   perimeterDiff = perimeter.uv - perimeter.nouv,
		   perimeterDiffSE = sqrt(perimeterSE.uv^2 + perimeterSE.nouv^2)
		   ) %>%
	select(mouse,Metadata_Time,eduDiff,caspDiff,tunelDiff,eduAreaDiff,caspAreaDiff,tunelAreaDiff,
		   eduDiffSE, caspDiffSE, tunelDiffSE, eduAreaDiffSE, caspAreaDiffSE, tunelAreaDiffSE, areaDiff, areaDiffSE, perimeterDiff, perimeterDiffSE)

diffGraph <- melt(diffs,id.vars=c("mouse","Metadata_Time"))
diffGraph$variable <- factor(diffGraph$variable, levels=c("eduDiff","caspDiff","tunelDiff","eduAreaDiff","caspAreaDiff","tunelAreaDiff","areaDiff","perimeterDiff"),
                                 labels=c("EdU","Caspase-3","TUNEL","EdU Area","Caspase-3 Area","TUNEL Area", "Area", "Perimeter"))


doTimePlots(diffGraph, "EdU",       "IR Intensity - NOIR Intensity", "EDU",       "eduIntensity.png")
doTimePlots(diffGraph, "Caspase-3", "IR Intensity - NOIR Intensity", "Caspase-3", "caspaseIntensity.png")
doTimePlots(diffGraph, "TUNEL",     "IR Intensity - NOIR Intensity", "TUNEL",     "tunelIntensity.png")

doTimePlots(diffGraph, "EdU Area",       "IR Area - NOIR Area", "EDU",       "eduArea.png")
doTimePlots(diffGraph, "Caspase-3 Area", "IR Area - NOIR Area", "Caspase-3", "caspaseArea.png")
doTimePlots(diffGraph, "TUNEL Area",     "IR Area - NOIR Area", "TUNEL",     "tunelArea.png")

doTimePlots(diffGraph, "Area",           "IR Area - NOIR Area", "Area",     "area.png")
doTimePlots(diffGraph, "Perimeter",     "IR Perimeter - NOIR Perimeter", "Perimeter",     "perimeter.png")

# Plots for each mouse
#controlDiffs <- filter(diffs, genotype == "Control")
#genoDiffs <- filter(diffs, genotype != "Control")

doMousePlots(diffs, "IR Intensity - NOIR Intensity", "Caspase-3 - Control", "caspDiff", "caspDiffSE", "caspaseIntensityAllMice-Control.png")
doMousePlots(diffs, "IR Intensity - NOIR Intensity", "EdU - Control", "eduDiff", "eduDiffSE", "eduIntensityAllMice-Control.png")
doMousePlots(diffs, "IR Intensity - NOIR Intensity", "TUNEL - Control", "tunelDiff", "tunelDiffSE", "tunelIntensityAllMice-Control.png")

doMousePlots(diffs, "IR Area - NOIR Area", "Caspase-3 - Control", "caspAreaDiff", "caspAreaDiffSE", "caspaseAreaAllMice-Control.png")
doMousePlots(diffs, "IR Area - NOIR Area", "EdU - Control", "eduAreaDiff", "eduAreaDiffSE", "eduAreaAllMice-Control.png")
doMousePlots(diffs, "IR Area - NOIR Area", "TUNEL - Control", "tunelAreaDiff", "tunelAreaDiffSE", "tunelAreaAllMice-Control.png")

doMousePlots(diffs, "IR Area - NOIR Area", "Area", "areaDiff", "areaDiffSE", "AreaAllMice-Control.png")
doMousePlots(diffs, "IR Perimeter - NOIR Perimeter", "Perimeter", "perimeterDiff", "perimeterDiffSE", "PerimeterAllMice-Control.png")

# Plot all the data
for (i in c("casp","caspArea","edu","eduArea","tunel","tunelArea")) {
	for (t in c("0HR","3HR","6HR","12HR","24HR")) {
		comparePlotsByTime(myData, i, t)
	}
}

# Do the anova's and write out the tukey's
# TODO: Change the label to make the csv readable
pVals <- data.frame("label"=c("Metadata_Time"))
writeHeader <- TRUE
for (x in c("eduDiff","caspDiff","tunelDiff","eduAreaDiff","caspAreaDiff","tunelAreaDiff")) {
  a <- aov(get(x) ~ Metadata_Time, data=diffs)
  pVals[x] <- summary(a)[[1]][["Pr(>F)"]][1]
  tky = cbind(x, TukeyHSD(a)$`Metadata_Time`)
  write.table(tky, file=paste(OUTDIR,"/tukey.csv", sep=""), append=!writeHeader, col.names=writeHeader, quote=FALSE, sep=",")
  writeHeader <- FALSE
}

# t-test for each variable and time point
out <- NULL
for (var in c("AreaShape_Area", "AreaShape_Perimeter", "CASP.intensity.corrected", "EDU.Intensity.corrected", "TUNEL.Intensity.Corrected", "CASP.area.corrected", "EDU.area.corrected", "TUNEL.area.Corrected")) {
	for (t in c("0HR","3HR","6HR","12HR","24HR")) {
		result <- t.test(formula = get(var) ~ Metadata_Treatment, data=dat, subset = Metadata_Time == t)
		out <- rbind(out, data.frame(var=var, t=t, pVal=result$p.value))
	}
}
write.table(out, file=paste(OUTDIR,"/pValuesByTime.csv", sep=""), append=FALSE, col.names=c("Measurement","Time","Pval"), row.names=FALSE, quote=FALSE, sep=",")

# t-test for each variable, time point, and mouse
out <- NULL
mouseTime <- distinct(dat, Metadata_Mouse, Metadata_Time)
for (var in c("AreaShape_Area", "AreaShape_Perimeter", "CASP.intensity.corrected", "EDU.Intensity.corrected", "TUNEL.Intensity.Corrected", "CASP.area.corrected", "EDU.area.corrected", "TUNEL.area.Corrected")) {
	for ( mT in 1:nrow(mouseTime)) {
		t <- mouseTime$Metadata_Time[mT]
		mouse <- mouseTime$Metadata_Mouse[mT]
		try(out <- rbind(out, data.frame(var=var, mouse=mouse, t=t,
			 pVal=t.test(formula = get(var) ~ Metadata_Treatment, data=dat, subset = Metadata_Time == t & Metadata_Mouse == mouse)$p.value)))
	}
}
write.table(out, file=paste(OUTDIR,"/pValuesByMouseAndTime.csv", sep=""), append=FALSE, col.names=c("Measurement","Mouse","Time","Pval"), row.names=FALSE, quote=FALSE, sep=",")



# Functions...

doTimePlots <- function(allDiffData,myVar,myYlabel,myTitle,fileName,myDodge=0.4) {

	diffGraph <- filter(allDiffData, variable == myVar)

    diffError <- diffGraph %>%
      group_by(Metadata_Time) %>%
      summarise( m = mean(value), se = sd(value)/sqrt(n()) )

    ggplot() + 
	  geom_hline(aes(yintercept=0), linetype="dotted") +
	  geom_line(aes(x = Metadata_Time, y = m), size = 1.5, linetype = 2, position=position_dodge(myDodge), data=diffError) +
      geom_dotplot(aes(y = value, x = Metadata_Time), data = diffGraph,
                   binaxis='y', stackdir='center', dotsize=1, stackgroups = TRUE, position=position_dodge(myDodge)) +
      geom_errorbar(aes(x = Metadata_Time, ymin=m-se, ymax=m), width=.2, position=position_dodge(myDodge), data=diffError) +
      geom_errorbar(aes(x = Metadata_Time, ymin=m, ymax=m+se), width=.2, position=position_dodge(myDodge), data=diffError) +
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
	  scale_x_discrete(limits=c("0HR","3HR","6HR","12HR","24HR")) +
      scale_fill_manual(values=myColors, name="", labels=myLabels) +
      scale_colour_manual(values=myColors, name="", labels=myLabels)
      ggsave(filename=fileName, path = OUTDIR, width=10, height=7, bg="transparent")
}


#allData <- myData
#colName <- "casp"
#tp <- "3HR"

comparePlotsByTime <- function(allData, colName, tp) {
 compData <- myData %>% filter(Metadata_Time == tp)
 fileName <- paste(colName,tp,"png",sep=".")
 ggplot(compData, aes(x = treatment, y = get(colName), label=mouse, group = mouse )) +
   geom_line(size=1.5) +
   geom_text(colour="black") +
   #geom_point(colour="black",pch=21, size=5) +
   #facet_wrap(~ genotype, labeller=label_parsed) +
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
	  geom_line(aes(x = Metadata_Time, y = get(colName), group=mouse, colour=mouse), size = 1.5, linetype = 2, position=position_dodge(myDodge), data=myData) +
      geom_errorbar(aes(x = Metadata_Time, fill=mouse, colour=mouse, ymin=get(colName)-get(errName), ymax=get(colName)), width=.2, position=position_dodge(myDodge), data=myData) +
      geom_errorbar(aes(x = Metadata_Time, fill=mouse, colour=mouse, ymin=get(colName), ymax=get(colName)+get(errName)), width=.2, position=position_dodge(myDodge), data=myData) +
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
	  scale_x_discrete(limits=c("0HR","3HR","6HR","12HR","24HR"))
      ggsave(filename=fileName, path = OUTDIR, width=10, height=7, bg="transparent")
}
