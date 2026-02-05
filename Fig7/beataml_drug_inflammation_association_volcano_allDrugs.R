library(data.table)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggthemes)

allDrug_Stat <- fread("./beataml_drug_inflammation_association_allDrugs.csv", data.table =F)

## volcano plot
plot(allDrug_Stat$LFC,-log2(allDrug_Stat$p))

library(EnhancedVolcano)
library(ggthemes)

pdf("beataml_drug_inflammation_association_volcano_allDrugs.pdf", width = 6, height = 5)
EnhancedVolcano(allDrug_Stat,
    lab = rownames(allDrug_Stat),
    x = 'LFC',
    y = 'p',
    xlab = "Log Fold Change (LFC)",
    ylab = bquote(~-Log[10] ~ italic(P)),
    pCutoff = 10e-6,
    FCcutoff = 0.5,
    pointSize = 2.5,
    labSize = 4,
    colAlpha = 0.7,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.75) +
    theme_few()
dev.off()
