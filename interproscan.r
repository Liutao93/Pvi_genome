#导入R包

library(tidyverse)

library(ggpubr)

# 读取结果文件并去掉有缺失值的行

interpro <- read_tsv("all.interpro",na = "N/A") %>% na.omit()

# 统计蛋白质家族、结构域和功能位点的比例等

ipr <- interpro %>% select(model,ipr_acc,ipr_desc) %>% group_by(model, ipr_acc) %>%
  
  summarise(ipr_desc = ipr_desc[[1]]) %>% group_by(ipr_acc, ipr_desc) %>% summarise(Count=n())%>%
  
  arrange(desc(Count)) %>% ungroup() %>%mutate(Percent = Count/sum(Count))

# 绘制前20的结构功能域或者功能位点

p <- ggplot(ipr) +
  
  geom_bar(aes(x = ipr_desc, y = Percent, fill = ipr_desc), stat = "identity") +
  
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.08),name = "Percent of Domain") +
  
  scale_x_discrete(limits = ipr$ipr_desc[1:20], name = NULL) + scale_fill_discrete(guide = FALSE)+
  
  theme_pubr() +
  
  theme(axis.text.x=element_text(angle=60,vjust=1, hjust=1))

#展示图片

p

#保存图片

ggsave("interpro.pdf", p, width = 16, height = 10)

ggsave("interpro.png", p, width = 16, height = 10)
