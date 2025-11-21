library(tidyverse)

pca <- read_table("/home/rymbekovaa95/Desktop/African/Afr1.eigenvec", col_names = FALSE)
eigenval <- scan("/home/rymbekovaa95/Desktop/African/Afr1.eigenval")

pca <- pca[,-1] #maybe no need

spp <- rep(NA, length(pca$X2))
spp

spp[grep("KSP062", pca$X2)] <- "Karretjie"
spp[grep("KSP063", pca$X2)] <- "Karretjie"
spp[grep("KSP065", pca$X2)] <- "Karretjie"
spp[grep("KSP067", pca$X2)] <- "Karretjie"
spp[grep("KSP069", pca$X2)] <- "Karretjie"

spp[grep("KSP092", pca$X2)] <- "Khutse"
spp[grep("KSP096", pca$X2)] <- "Khutse"

spp[grep("KSP103", pca$X2)] <- "Ju'hoansi"
spp[grep("KSP105", pca$X2)] <- "Ju'hoansi"
spp[grep("KSP106", pca$X2)] <- "Ju'hoansi"
spp[grep("KSP111", pca$X2)] <- "Ju'hoansi"
spp[grep("KSP116", pca$X2)] <- "Ju'hoansi"

spp[grep("KSP124", pca$X2)] <- "Nama"
spp[grep("KSP134", pca$X2)] <- "Nama"
spp[grep("KSP137", pca$X2)] <- "Nama"
spp[grep("KSP139", pca$X2)] <- "Nama"
spp[grep("KSP140", pca$X2)] <- "Nama"

spp[grep("KSP146", pca$X2)] <- "Xun"
spp[grep("KSP150", pca$X2)] <- "Xun"
spp[grep("KSP152", pca$X2)] <- "Xun"
spp[grep("KSP154", pca$X2)] <- "Xun"
spp[grep("KSP155", pca$X2)] <- "Xun"

spp[grep("KSP224", pca$X2)] <- "Khutse San"
spp[grep("KSP225", pca$X2)] <- "Khutse San"
spp[grep("KSP228", pca$X2)] <- "Khutse San"


spp[grep("PNP000", pca$X2)] <- "Baka"
spp[grep("PNP001", pca$X2)] <- "Baka"
spp[grep("PNP002", pca$X2)] <- "Baka"
spp[grep("PNP003", pca$X2)] <- "Baka"
spp[grep("PNP004", pca$X2)] <- "Baka"
spp[grep("PNP005", pca$X2)] <- "Baka"
spp[grep("PNP006", pca$X2)] <- "Baka"

spp[grep("PNP010", pca$X2)] <- "Nzime"
spp[grep("PNP011", pca$X2)] <- "Nzime"
spp[grep("PNP012", pca$X2)] <- "Nzime"
spp[grep("PNP013", pca$X2)] <- "Nzime"
spp[grep("PNP014", pca$X2)] <- "Nzime"

spp[grep("PNP020", pca$X2)] <- "Ba.Kola"
spp[grep("PNP021", pca$X2)] <- "Ba.Kola"
spp[grep("PNP022", pca$X2)] <- "Ba.Kola"
spp[grep("PNP023", pca$X2)] <- "Ba.Kola"
spp[grep("PNP024", pca$X2)] <- "Ba.Kola"

spp[grep("PNP030", pca$X2)] <- "Ngumba"
spp[grep("PNP031", pca$X2)] <- "Ngumba"
spp[grep("PNP032", pca$X2)] <- "Ngumba"
spp[grep("PNP033", pca$X2)] <- "Ngumba"
spp[grep("PNP034", pca$X2)] <- "Ngumba"

spp[grep("PNP040", pca$X2)] <- "Aka Mbati"
spp[grep("PNP041", pca$X2)] <- "Aka Mbati"
spp[grep("PNP042", pca$X2)] <- "Aka Mbati"
spp[grep("PNP043", pca$X2)] <- "Aka Mbati"
spp[grep("PNP044", pca$X2)] <- "Aka Mbati"

spp[grep("PNP050", pca$X2)] <- "Ba.Kiga"
spp[grep("PNP051", pca$X2)] <- "Ba.Kiga"
spp[grep("PNP052", pca$X2)] <- "Ba.Kiga"
spp[grep("PNP053", pca$X2)] <- "Ba.Kiga"
spp[grep("PNP054", pca$X2)] <- "Ba.Kiga"

spp[grep("PNP060", pca$X2)] <- "Ba.Twa"
spp[grep("PNP062", pca$X2)] <- "Ba.Twa"
spp[grep("PNP063", pca$X2)] <- "Ba.Twa"
spp[grep("PNP064", pca$X2)] <- "Ba.Twa"
spp[grep("PNP065", pca$X2)] <- "Ba.Twa"
spp[grep("PNP066", pca$X2)] <- "Ba.Twa"

spp[grep("PNP070", pca$X2)] <- "Nsua"
spp[grep("PNP071", pca$X2)] <- "Nsua"
spp[grep("PNP072", pca$X2)] <- "Nsua"
spp[grep("PNP073", pca$X2)] <- "Nsua"
spp[grep("PNP074", pca$X2)] <- "Nsua"

spp[grep("PNP080", pca$X2)] <- "Ba.Konjo"
spp[grep("PNP081", pca$X2)] <- "Ba.Konjo"
spp[grep("PNP082", pca$X2)] <- "Ba.Konjo"
spp[grep("PNP083", pca$X2)] <- "Ba.Konjo"
spp[grep("PNP084", pca$X2)] <- "Ba.Konjo"

pca <- as_tibble(data.frame(pca, spp))
pca <- pca[complete.cases(pca), ]
spp <- spp[complete.cases(spp)]

spp1 <- rep(NA, length(73))
spp1

spp1 <- ifelse(pca$X2 %in% c("KSP062", "KSP063", "KSP065", "KSP067", "KSP069", "KSP092", "KSP096", "KSP103", "KSP105", "KSP106", "KSP111", "KSP116", "KSP124", "KSP134", "KSP137", "KSP139", "KSP140", "KSP146", "KSP150", "KSP152", "KSP154", "KSP155", "KSP224", "KSP225", "KSP228"), "Khoesan (KS)", spp1)

spp1 <- ifelse(pca$X2 %in% c("PNP000", "PNP001", "PNP002", "PNP003", "PNP004", "PNP005", "PNP006", "PNP020", "PNP021", "PNP022", "PNP023", "PNP024", "PNP040", "PNP041", "PNP042", "PNP043", "PNP044", "PNP060", "PNP062", "PNP063", "PNP064", "PNP065", "PNP066", "PNP070", "PNP071", "PNP072", "PNP073", "PNP074"), "rainforest hunter-gatherers (RHG)", spp1)

spp1 <- ifelse(pca$X2 %in% c("PNP010", "PNP011", "PNP012", "PNP013", "PNP014", "PNP030", "PNP031", "PNP032", "PNP033", "PNP034", "PNP050", "PNP051", "PNP052", "PNP053", "PNP054", "PNP080", "PNP081", "PNP082", "PNP083", "PNP084"), "rainforest hunter-gatherers neighbors (RHGn)", spp1)



spp1

pca <- as_tibble(data.frame(pca, spp, spp1))

pca


pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)


b <- ggplot(pca, aes(X3, X4, col = spp, shape = spp1)) + geom_point(size = 5)
b <- b + coord_equal() + theme_light()
b <- b + labs(col = "Population", shape = "Sample")
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b
#variance explained
b <- ggplot(pca, aes(X3, X4, col = spp, shape = spp1)) + geom_point(size = 5)
b <- b + coord_equal() + theme_light()
b <- b + labs(col = "Population", shape = "Sample")

# Updated labels for the axes
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "% Variance Explained)")) 
b <- b + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "% Variance Explained)"))

# Display the plot
b


ggsave("African_names_colors.png", plot = b, width = 10, height = 8, dpi = 300)
