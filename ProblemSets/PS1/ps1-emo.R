library(ggplot2)

df = read.csv('/Users/erin/Desktop/GitHub/emochoa/MACS30200proj/ProblemSets/PS1/df.csv')
#dfal = read.csv('/Users/erin/Desktop/GitHub/emochoa/MACS30200proj/ProblemSets/PS1/dfal.csv')
#dfbh = read.csv('/Users/erin/Desktop/GitHub/emochoa/MACS30200proj/ProblemSets/PS1/dfbh.csv')

df$EduCat[df$MaternalEducation == '8th grade or less'] = 'Associate\'s or lower'
df$EduCat[df$MaternalEducation == 'Associate\'s'] = 'Associate\'s or lower'
df$EduCat[df$MaternalEducation == 'High school, no diploma'] = 'Associate\'s or lower'
df$EduCat[df$MaternalEducation == 'HSD/GED'] = 'Associate\'s or lower'
df$EduCat[df$MaternalEducation == 'Some college, no degree'] = 'Associate\'s or lower'
df$EduCat[df$MaternalEducation == 'Bachelor\'s'] = 'Bachelor\'s or higher'
df$EduCat[df$MaternalEducation == 'Doctorate/Professional'] = 'Bachelor\'s or higher'
df$EduCat[df$MaternalEducation == 'Master\'s'] = 'Bachelor\'s or higher'

df = df[(!is.na(df$EduCat)), ]

df$EduCat = as.factor(df$EduCat)

ggplot(df, mapping = aes(x = MaternalAge, fill = EduCat)) +
  scale_fill_discrete(name="Education Level") + 
  geom_histogram(alpha = .7, binwidth = 2) +
  geom_vline(xintercept=almean, color = 'coral') + geom_vline(xintercept=bhmean, color='darkturquoise') +
  labs(title = "Maternal Age by Maternal Education Level",
       x = "Maternal age",
       y = "Frequency count of birth mothers") +
  theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(linetype = "solid", color = "grey70", fill=NA, size=1.1)) +
  scale_x_continuous(breaks = seq(12,52,by=4))


