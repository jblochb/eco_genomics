library(tidyverse)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

Met<- read_csv("ThermalStressMetabolic.csv")
setwd ("/gpfs1/cl/pbio3990/Intro_to_R/")
Met
dim(Met)
head(Met)
str(Met)
colnames(Met)
view(Met)
Met %>%
  rename(SMR = 'SMR (g/hr)')  
Met_trans <- Met %>%
  rename(SMR = 'SMR (g/hr)', MMR = 'MMR (g/hr)')
Met_trans2 <- Met_trans
  rename(MMR = 'MMR (g/hr)')
view(Met_trans)

Met_trans %>%
  select(Temp, Line, SMR, MMR)
Met_trans %>%
  mutate(ThermalScope = MMR -SMR)
colnames(Met_trans)
Met_trans %>%
  filter(Temp == 34 & Line == "LS1")
view(Met_trans %>%
  filter(Temp ==34 | Temp == 28))
#arrange data from lowest to highest or highest to lowest 
Met_trans %>%
  arrange(SMR)
Met_trans %>%
  arrange(desc(SMR))
#to identify any duplicates
Met_trans %>%
  distinct()
Met_trans %>%
  rename(Temp.celcius = "Temp")
Met_new <- Met_trans %>%
  arrange(desc(r2)) %>%
  mutate(ThermalScope = MMR - SMR)
#To see which the highest r2 value for each temperature 
view(Met_new %>%
       distinct(Temp, r2))
Met_new %>%
  group_by(Temp) %>%
  summarize(avg_r2 = mean(r2, na.rm = T))
#to see a few more state
Met_new %>%
  group_by(Temp) %>%
  drop_na() %>%
  summarise(
    mean=mean(SMR),
    sd=sd(SMR),
    min=min(SMR),
    max=max(SMR),
    n=n())
#see which treatment has the highest MMR
Met_new %>%
  arrange(desc(MMR)) %>% select(MMR, Temp, Treat, Line)
#count different variations in color
view(Met_new %>% distinct(Color))

#plotting

library(ggplot2)
ggplot(Met_trans, aes(x = Temp, y = SMR, fill = Temp))+
  geom_boxplot()+
  scale_fill_manual(values = c("blue", "orange", "red"))
ggplot(Met_trans, aes(x = Temp, y = SMR, fill = Temp))+
  geom_boxplot()+
  scale_fill_brewer(palette = "Dark2")
Met_trans$Temp <- as.factor(Met_trans$Temp)
ggplot(Met_trans, aes(x=MMR, y=SMR, colour=Temp))+
  geom_point()+
  scale_color_colorblind()+
  facet_wrap(~Line)+
  labs(title = "Practice R Graph", x= "MMR (units)")+
  theme(text=element_text(size = 16, family = "TT Comic Sans MS"))
library(ggthemes)
install.packages("ggthemes")
