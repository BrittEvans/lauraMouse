library(dplyr)
library(ggplot2)
library(reshape2)

dat <- read.csv("~/lauraMouse/data/tg.csv")
dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c(expression("Villin-CRE BAI1"^{WT/WT}),expression("Villin-CRE BAI1"^{WT/Tg}))
myLevels <- c("WT","Tg")
myColors <- c("yellow","cyan")
OUTDIR <- "~/lauraMouse/tg"
doitAll()

dat <- read.csv("~/lauraMouse/data/belmo.csv")
myLabels <- c(expression("Villin-CRE BELMO"^{WT/WT}),expression("Villin-CRE BELMO"^{WT/Tg}))
myLevels <- c("WT","BELMO")
myColors <- c("yellow","green")
OUTDIR <- "~/lauraMouse/belmo"
doitAll()

dat <- read.csv("~/lauraMouse/data/cko.csv")
myLabels <- c(expression("Villin-CRE BAI1"^{FL/FL}),expression("Villin-CRE BAI1"^{Delta/Delta}))
myLevels <- c("WT","KO")
myColors <- c("yellow","red")
OUTDIR <- "~/lauraMouse/cko"
doitAll()

dat <- read.csv("~/lauraMouse/data/bai1.csv")
dat <- dat %>% rename(Metadata_Frame = GENOTYPE, Metadata_Treatment = TREATMENT, Metadata_Mouse = MOUSE..)
dat <- dat %>% filter(Metadata_Mouse != 'AV3')
myLabels <- c(expression("Villin-CRE BAI1"^{FL/FL}),expression("Villin-CRE BAI1"^{Delta/Delta}))
myLevels <- c("WT","KO")
myColors <- c("yellow","red")
OUTDIR <- "~/lauraMouse/bai1"
doitAll()

dat <- read.csv("~/lauraMouse/data/tg2.csv")
dat <- dat %>% rename(Metadata_Frame = GENOTYPE, Metadata_Treatment = TREATMENT, Metadata_Mouse = MOUSE..)
myLabels <- c(expression("Villin-CRE BAI1"^{WT/WT}),expression("Villin-CRE BAI1"^{WT/Tg}))
myLevels <- c("WT","Tg")
myColors <- c("yellow","cyan")
OUTDIR <- "~/lauraMouse/tg2"
doitAll()


doitAll <- function() {
  # Read and group the data
  grouped <- dat %>%
    #mutate(genotype = factor(Metadata_Genotype, levels = myLevels)) %>%
    mutate(genotype = factor(Metadata_Frame, levels = myLevels, labels=myLabels)) %>%
    rename(mouse = Metadata_Mouse, treatment = Metadata_Treatment) %>%
    group_by(mouse,treatment,genotype) %>%
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
               tunelSE = sd(TUNEL.Intensity.Corrected)/sqrt(count) )
  
  # Write out the csv
  write.csv(grouped, paste(OUTDIR,"metrics.csv",sep="/"))
  
  # Do the anova's
  tukeys <- data.frame("label"=c("variable","diff","lwr","upr","pVal"))
  pVals <- data.frame("label"=c("genotype","treatment","interaction"))
  writeHeader <- TRUE
  for (x in c("area","perimeter","casp","edu","tunel","count")) {
    a <- aov(get(x) ~ genotype*treatment, data=grouped)
    pVals[x] <- summary(a)[[1]][["Pr(>F)"]][1:3]
    tky = cbind(x, TukeyHSD(a)$`genotype:treatment`)
    write.table(tky, file=paste(OUTDIR,"/tukey.csv", sep=""), append=!writeHeader, col.names=writeHeader, quote=FALSE, sep=",")
    writeHeader <- FALSE
  }
  write.csv(pVals,paste(OUTDIR,"anova.csv",sep="/"))
  
  makeDotplots( grouped %>% select (casp, treatment, genotype) %>% rename(myColumn=casp), "Caspase-3 Intensity/hoechst Intensity", "casp.png")
  newPlotsWithBars( grouped %>% select (casp, treatment, genotype) %>% rename(myColumn=casp), "Caspase-3 Intensity/hoechst Intensity", "casp-panels-bars.png")
  newPlotsWithoutBars( grouped %>% select (casp, treatment, genotype) %>% rename(myColumn=casp), "Caspase-3 Intensity/hoechst Intensity", "casp-panels-noBars.png")

  makeDotplots( grouped %>% select (edu, treatment, genotype) %>% rename(myColumn=edu), "EdU Intensity/hoechst Intensity", "edu.png")
  newPlotsWithBars( grouped %>% select (edu, treatment, genotype) %>% rename(myColumn=edu), "EdU Intensity/hoechst Intensity", "edu-panels-bars.png")
  newPlotsWithoutBars( grouped %>% select (edu, treatment, genotype) %>% rename(myColumn=edu), "EdU Intensity/hoechst Intensity", "edu-panels-noBars.png")

  makeDotplots( grouped %>% select (tunel, treatment, genotype) %>% rename(myColumn=tunel), "TUNEL Intensity/hoechst Intensity", "tunel.png")
  newPlotsWithBars( grouped %>% select (tunel, treatment, genotype) %>% rename(myColumn=tunel), "TUNEL Intensity/hoechst Intensity", "tunel-panels-bars.png")
  newPlotsWithoutBars( grouped %>% select (tunel, treatment, genotype) %>% rename(myColumn=tunel), "TUNEL Intensity/hoechst Intensity", "tunel-panels-noBars.png")

  makeDotplots( grouped %>% select (area, treatment, genotype) %>% rename(myColumn=area), "Area", "area.png")
  newPlotsWithBars( grouped %>% select (area, treatment, genotype) %>% rename(myColumn=area), "Area", "area-panels-bars.png")
  newPlotsWithoutBars( grouped %>% select (area, treatment, genotype) %>% rename(myColumn=area), "Area", "area-panels-noBars.png")

  makeDotplots( grouped %>% select (perimeter, treatment, genotype) %>% rename(myColumn=perimeter), "Perimeter", "perimeter.png")
  newPlotsWithBars( grouped %>% select (perimeter, treatment, genotype) %>% rename(myColumn=perimeter), "Perimeter", "perimeter-panels-bars.png")
  newPlotsWithoutBars( grouped %>% select (perimeter, treatment, genotype) %>% rename(myColumn=perimeter), "Perimeter", "perimeter-panels-noBars.png")

  makeDotplots( grouped %>% select (count, treatment, genotype) %>% rename(myColumn=count), "Count", "count.png")
  newPlotsWithBars( grouped %>% select (count, treatment, genotype) %>% rename(myColumn=count), "Count", "count-panels-bars.png")
  newPlotsWithoutBars( grouped %>% select (count, treatment, genotype) %>% rename(myColumn=count), "Count", "count-panels-noBars.png")
  
  diffPlots( grouped, "intensityDiff.png")
}


# Make graphs
makeDotplots <- function(myData, yLabel, fileName) {

  grouped2 <- myData %>%
    group_by(genotype,treatment) %>%
    summarise( m = mean(myColumn), se = sd(myColumn)/sqrt(n()) )

  ggplot() + 
    geom_dotplot(aes(y = myColumn, x = treatment, fill = genotype), data = myData,
                 binaxis='y', stackdir='center', dotsize=1, stackgroups = TRUE, position="dodge") +
    geom_errorbar(aes(x = treatment,  fill = genotype, ymin=m-se, ymax=m), width=.2, position=position_dodge(.9), data=grouped2) +
    geom_errorbar(aes(x = treatment,  fill = genotype, y=m, ymin=m, ymax=m+se), width=.2, position=position_dodge(.9), data=grouped2) +
    theme_classic(base_size=14) +
    theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.text.align = 0,
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          rect = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", colour= NA),
          panel.background = element_rect(fill = "transparent", colour= NA),
          legend.text=element_text(size=14)) +
    labs(y=yLabel,x="Treatment") +
    scale_fill_manual(name="Genotype", labels=myLabels, values=myColors)
  ggsave(paste(OUTDIR, fileName, sep="/"), bg="transparent")
}


# Differences
#UV - NOUV per mouse
#casp edu tunel
# y-axis: UV Intensity - NOUV Intensity
diffPlots <- function(myData, fileName) {
    uv <- filter(myData, treatment == "UV")
    nouv <- filter(myData, treatment == "NOUV")
    diffs <- inner_join(uv, nouv, by=c("mouse","genotype"), suffix=c(".uv",".nouv")) %>%
        ungroup() %>%
        mutate(eduDiff = edu.uv - edu.nouv,
               caspDiff = casp.uv - casp.nouv,
               tunelDiff = tunel.uv - tunel.nouv) %>%
        select(mouse,genotype,eduDiff,caspDiff,tunelDiff)
  
    diffGraph <- melt(diffs)
    diffGraph$variable <- factor(diffGraph$variable, levels=c("eduDiff","caspDiff","tunelDiff"),
                                 labels=c("EdU","Caspase-3","TUNEL"))
    diffError <- diffGraph %>%
      group_by(genotype,variable) %>%
      summarise( m = mean(value), se = sd(value)/sqrt(n()) )

    ggplot() + 
      geom_dotplot(aes(y = value, x = variable, fill = genotype), data = diffGraph,
                   binaxis='y', stackdir='center', dotsize=1, stackgroups = TRUE, position="dodge") +
      geom_errorbar(aes(x = variable,  fill = genotype, ymin=m-se, ymax=m), width=.2, position=position_dodge(.9), data=diffError) +
      geom_errorbar(aes(x = variable,  fill = genotype, y=m, ymin=m, ymax=m+se), width=.2, position=position_dodge(.9), data=diffError) +
      theme_classic(base_size=14) +
      theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
            legend.text.align = 0,
            axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
            rect = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent", colour= NA),
            panel.background = element_rect(fill = "transparent", colour= NA),
            legend.text=element_text(size=14)) +
      labs(y="UV Intensity - NOUV Intensity",x="") +
      scale_fill_manual(values=myColors, name="Genotype", labels=myLabels)
      ggsave(filename=fileName, path = OUTDIR, width=10, height=7, bg="transparent")

    output <- data.frame("variable" = c("edu","casp","tunel"))
    output["pValue"] <- c(t.test(diffs$eduDiff)$p.value, t.test(diffs$caspDiff)$p.value, t.test(diffs$tunelDiff)$p.value)
    write.csv(output,paste(OUTDIR,"intensityDifference.csv",sep="/"))
}

newPlotsWithBars <- function(myData, yLabel, fileName) {

  grouped2 <- myData %>%
    group_by(genotype,treatment) %>%
    summarise( m = mean(myColumn), se = sd(myColumn)/sqrt(n()) )

 ggplot() +
   geom_line(aes(x=treatment, y=myColumn, colour=genotype, group=mouse), data=myData, size=1.5) +
   geom_point(aes(x=treatment, y=myColumn, colour=genotype, fill=genotype), colour="black",pch=21, size=5, data=myData) +
    geom_errorbar(aes(x=treatment, ymin=m-se, ymax=m), width=.2, data=grouped2) +
    geom_errorbar(aes(x=treatment, ymin=m, ymax=m+se), width=.2, data=grouped2) +
   facet_wrap(~ genotype, labeller=label_parsed) +
   scale_colour_manual(name="Genotype", values=myColors, guide=FALSE) +
   scale_fill_manual(name="Genotype", values=myColors, guide=FALSE) +
    theme_classic(base_size=14) +
    theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.text.align = 0,
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          rect = element_rect(fill = "transparent"),
		  strip.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour= NA),
          panel.background = element_rect(fill = "transparent", colour= NA),
          legend.text=element_text(size=14)) +
    labs(y=yLabel,x="Treatment")
  ggsave(paste(OUTDIR, fileName, sep="/"), bg="transparent")
}

newPlotsWithoutBars <- function(myData, yLabel, fileName) {

 ggplot() +
   geom_line(aes(x=treatment, y=myColumn, colour=genotype, group=mouse), data=myData, size=1.5) +
   geom_point(aes(x=treatment, y=myColumn, colour=genotype, fill=genotype), colour="black",pch=21, size=5, data=myData) +
   facet_wrap(~ genotype, labeller=label_parsed) +
   scale_colour_manual(name="Genotype", values=myColors, guide=FALSE) +
   scale_fill_manual(name="Genotype", values=myColors, guide=FALSE) +
    theme_classic(base_size=14) +
    theme( axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          legend.text.align = 0,
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          rect = element_rect(fill = "transparent"),
		  strip.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour= NA),
          panel.background = element_rect(fill = "transparent", colour= NA),
          legend.text=element_text(size=14)) +
    labs(y=yLabel,x="Treatment")
  ggsave(paste(OUTDIR, fileName, sep="/"), bg="transparent")
}
