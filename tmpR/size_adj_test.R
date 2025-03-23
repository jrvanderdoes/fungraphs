results_size <- readRDS(paste0('C:/Users/jerem/OneDrive/Documents/',
                         'School/Waterloo/Research/',
                         'GraphCP/Data/Simulations/power/',
                         'trials',trials,
                         '_ns',n,
                         '_bl',break_location,
                         '_dep_size',
                         '_results_simp.rds'))
results_size$mean <- 1-as.numeric(results_size$mean)
results_size$cov <- 1-as.numeric(results_size$cov)
data_size <- results_size %>% pivot_longer(cols = MDT1ori:cov)
data_size$value <- as.numeric(data_size$value)

results_power <- readRDS(paste0('C:/Users/jerem/OneDrive/Documents/',
                               'School/Waterloo/Research/',
                               'GraphCP/Data/Simulations/power/',
                               'trials',trials,
                               '_ns',n,
                               '_bl',break_location,
                               '_dep',
                               '_results_simp.rds'))
results_power$mean <- 1-as.numeric(results_power$mean)
results_power$cov <- 1-as.numeric(results_power$cov)
data_power <- results_power %>% pivot_longer(cols = MDT1ori:cov)
data_power$value <- as.numeric(data_power$value)

## Full_data
data_power1 <- data_power
colnames(data_power1) <- c('dep','name','power')
data_size1 <- data_size
colnames(data_size1) <- c('dep','name','size')

data_full <- merge(data_power1,data_size1,by = c('name','dep'))
data_full$sizeTo5 <- data_full$size-0.05
data_full$sizeAdjPower <- data_full$power-data_full$sizeTo5

########################
tmp <- ggplot() +
  geom_line(aes(x=dep,
                y=value,
                color=name,group=name),data_power) +
  geom_line(aes(x=dep,
                y=value,
                color=name,group=name),data_size,linetype='dashed') +
  geom_point(aes(x=dep,
                 y=value,
                 color=name,group=name),data_power) +
  geom_point(aes(x=dep,
                 y=value,
                 color=name,group=name),data_size) +
  geom_hline(aes(yintercept=0.05)) +
  geom_line(aes(x=dep,
                y=value,group=name),
            data=na.omit(data_power[data_power$name=='mean',]),
            color='black', linewidth=1,linetype='dashed') +
  geom_line(aes(x=dep,
                y=value,group=name),
            data=na.omit(data_size[data_size$name=='mean',]),
            color='black', linewidth=1,linetype='dashed') +
  geom_line(aes(x=dep,
                y=value,group=name),
            data=na.omit(data_power[data_power$name=='cov',]),
            color='black', linewidth=1, linetype='dotted') +
  geom_line(aes(x=dep,
                y=value,group=name),
            data=na.omit(data_size[data_size$name=='cov',]),
            color='black', linewidth=1, linetype='dotted') +
  guides(color="none") +
  theme_bw() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=18))+
  xlab('Strength of Dependence (mean change 0.25)') +
  ylab('Null Rejection Rate') +
  ylim(c(0,1)) #+
#xlim(10,1) #+
#scale_x_continuous(breaks=c(10:1),limits = c(10,1))

png("C:/Users/jerem/Downloads/dep_power_mean.png")
print(tmp)
dev.off()

########################

tmp <- ggplot(na.omit(data_full)) +
  geom_line(aes(x=dep,
                y=sizeAdjPower,
                color=name,group=name)) +
  geom_point(aes(x=dep,
                 y=sizeAdjPower,
                 color=name,group=name)) +
  geom_hline(aes(yintercept=0.05)) +
  geom_line(aes(x=dep,
                y=sizeAdjPower,group=name),
            data=na.omit(data_full[data_full$name=='mean',]),
            color='black', linewidth=1) +
  geom_line(aes(x=dep,
                y=sizeAdjPower,group=name),
            data=na.omit(data_full[data_full$name=='cov',]),
            color='black', linewidth=1, linetype='dashed') +
  guides(color="none") +
  theme_bw() +
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=18))+
  xlab('Strength of Dependence (mean change 0.25)') +
  ylab('Size Adjusted Rejection Rate') +
  ylim(c(0,1)) #+
#xlim(10,1) #+
#scale_x_continuous(breaks=c(10:1),limits = c(10,1))

png("C:/Users/jerem/Downloads/dep_size_adj_power.png")
print(tmp)
dev.off()

