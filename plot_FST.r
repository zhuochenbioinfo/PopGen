library(ggplot2)

args <- commandArgs(TRUE)
file <- args[1]
outfile <- args[2]

data <- read.delim(file)

goodChrOrder <- paste("Chr", c(1:12), sep="")
data$CHROM <- factor(data$CHROM, levels=goodChrOrder)

# 获取显著性阈值
top5.value <- as.numeric(quantile(data$WEIGHTED_FST,prob=1-5/100))
top10.value <- as.numeric(quantile(data$WEIGHTED_FST,prob=1-10/100))

# 根据显著性阈值加标签
data$sig <- "nosig"
data[data$WEIGHTED_FST >= top10.value,]$sig <- "top10pct"
data[data$WEIGHTED_FST >= top5.value,]$sig <- "top5pct"

# 将标签列转换为factor属性
data$sig <- as.factor(data$sig)

pdf(outfile,width=20,height=4)
ggplot(data[data$WEIGHTED_FST > 0,], aes(x=BIN_START/1000000, y=WEIGHTED_FST, color=sig)) + 
  geom_point(stat="identity", position="identity",size=0.5) + 
  # 设定颜色
  scale_color_manual(values = c("grey", "orange", "red")) +
  facet_wrap(~ CHROM, ncol=12, scales = "free_x", strip.position = "bottom") + 
  ggtitle("REF:MSU") +
  xlab("Position(Mb)") +
  ylab("value") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()
