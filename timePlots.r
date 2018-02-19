library(dplyr)
library(ggplot2)
library(reshape2)

dat <- read.csv("~/lauraMouse/timeplots/allTPdata.csv")
#dat <- dat %>% rename(Metadata_Frame = Metadata_Genotype)
myLabels <- c("Control",expression(italic("Villin-cre, Bai1"^{"WT/Tg"})))
myLevels <- c("WT","Tg")
myColors <- c("yellow","cyan")
OUTDIR <- "~/lauraMouse/timeplots/"

myData <- dat %>%
  mutate(genotype = factor(GENOTYPE, levels = myLevels, labels=myLabels)) %>%
  rename(mouse = MOUSE.., treatment = TREATMENT) %>%
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



# GENOTYPE (WT/Tg), MOUSE #, TIMEPOINT (3/6/12/24HR), TREATMENT (Un/nouv)

# 3 Total intensity
# 3 Total area corrected
# Repeat for CASPASE, TUNEL

# 2-way anova using time and variable
# And the tukey test



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
