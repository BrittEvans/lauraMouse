myPath <- c("abca1", "abca7", "adgrb1", "ager", "aif1", "alox15", "ano6", "appl2", "arhgap12", "arhgap25", "becn1", "bin2", "c3", "cd300a", "cd36", "clcn3", "clec7a", "elmo1", "f2rl1", "fcer1g", "fcgr1", "fcgr2b", "fcgr3", "gata2", "gsn", "gulp1", "igh-8", "igha", "ighd", "ighe", "ighg", "ighg1", "ighg2a", "ighg2b", "ighg2c", "ighg3", "ighm", "ighv1-11", "ighv1-12", "ighv1-15", "ighv1-16", "ighv1-18", "ighv1-19", "ighv1-20", "ighv1-22", "ighv1-23", "ighv1-24", "ighv1-26", "ighv1-31", "ighv1-34", "ighv1-36", "ighv1-37", "ighv1-39", "ighv1-4", "ighv1-42", "ighv1-43", "ighv1-47", "ighv1-49", "ighv1-5", "ighv1-50", "ighv1-52", "ighv1-53", "ighv1-54", "ighv1-55", "ighv1-56", "ighv1-58", "ighv1-59", "ighv1-61", "ighv1-62-1", "ighv1-62-2", "ighv1-62-3", "ighv1-63", "ighv1-64", "ighv1-66", "ighv1-67", "ighv1-69", "ighv1-7", "ighv1-72", "ighv1-74", "ighv1-75", "ighv1-76", "ighv1-77", "ighv1-78", "ighv1-80", "ighv1-81", "ighv1-82", "ighv1-84", "ighv1-85", "ighv1-9", "ighv10-1", "ighv10-3", "ighv11-1", "ighv11-2", "ighv12-3", "ighv13-2", "ighv14-1", "ighv14-2", "ighv14-3", "ighv14-4", "ighv15-2", "ighv16-1", "ighv2-2", "ighv2-3", "ighv2-4", "ighv2-5", "ighv2-6", "ighv2-6-8", "ighv2-7", "ighv2-9", "ighv2-9-1", "ighv3-1", "ighv3-3", "ighv3-4", "ighv3-5", "ighv3-6", "ighv3-8", "ighv4-1", "ighv5-12", "ighv5-12-4", "ighv5-15", "ighv5-16", "ighv5-17", "ighv5-2", "ighv5-4", "ighv5-6", "ighv5-9", "ighv5-9-1", "ighv6-3", "ighv6-4", "ighv6-5", "ighv6-6", "ighv6-7", "ighv7-2", "ighv7-4", "ighv8-11", "ighv8-12", "ighv8-13", "ighv8-2", "ighv8-4", "ighv8-5", "ighv8-6", "ighv8-8", "ighv8-9", "ighv9-1", "ighv9-2", "ighv9-3", "ighv9-4", "igkc", "iglc1", "iglc2", "iglc3", "igll1", "itga2", "itgb2", "lbp", "marco", "megf10", "mfge8", "msr1", "myh9", "nckap1l", "pparg", "rab31", "rac1", "rac3", "rhobtb1", "rhobtb2", "sh3bp1", "siglece", "sirpa", "stap1", "thbs1", "trbc1", "trbc2", "trdc", "trdv4", "trem2", "treml4", "vamp7", "xkr4", "xkr6", "xkr7", "xkr8", "xkr9")

dat <- read.csv("~/lauraMouse/rnaSeq/snrData.csv", header=TRUE)

head(dat)

myGenes <- dat %>%
	mutate(geneLow = tolower(gene)) %>%
	filter(geneLow %in% myPath)

write.table(myGenes, file = "~/lauraMouse/rnaSeq/filteredData.csv", sep=",", row.names=FALSE)

toPlot <- melt(myGenes, id.vars=c("gene")) %>%
	filter(variable %in% c("normSNR2","normSNR3")) %>%
	filter(value != "#DIV/0!")

#toPlot <- melt(myGenes, id.vars=c("gene")) %>%
#	filter(variable %in% c("snr2","snr3")) %>%
#	filter(value != "#DIV/0!")

ggplot(toPlot, aes(x=variable, y=gene)) +
	geom_tile(aes(fill=as.numeric(value))) +
	scale_fill_gradient(low = "red", high = "green")

	geom_point() +
	scale_colour_continuous()


