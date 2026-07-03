#####################################################
## summaries and maps of pesticide risk
## Authors: Francesca Mancini & Susan Jarvis
## Date created: 2025-06-03
#####################################################
options(scipen = 999)

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(terra)
library(tidyterra)
library(ggpubr)
library(scales)
library(rcartocolor)
library(patchwork)

ai_risk_overall <- readRDS("ai_risk_overall_corrected.rds")
ai_risk_by_chemgroup <- readRDS("ai_risk_chemgroup_corrected.rds")
ai_risk_by_type <- readRDS("ai_risk_type_corrected.rds")



# create maps of all raster stacks

all_files <- list.files("/data/outputs", pattern = ".tif")

for(f in 1:length(all_files)) {
  map <- rast(file.path("/data/outputs", all_files[f]))
  
  maps <- ggplot() +
    geom_spatraster(data = map) +
    facet_wrap(~lyr, ncol = 4) +
    scale_fill_viridis_c(na.value = NA,
                         name = "Risk") +
    ggtitle(gsub("_", " ", 
                 gsub(".tif", "", all_files[f]))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
  ggsave(plot = maps, width = 15, height = 20, units = "cm", dpi = 300,
         filename = file.path("/data/outputs",
                              paste0(gsub(".tif", "", all_files[f]), ".png")))
  
}


## How has the overall risk changed between 1994 and 2016 for each taxon?

ai_risk_summary <- ai_risk_by_chemgroup %>%
  group_by(year, chemgroup) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE),
            bumblebees_risk = sum(bumblebees_risk, na.rm = TRUE),
            obicornis_risk = sum(obicornis_risk, na.rm = TRUE),
            ocornuta_risk = sum(ocornuta_risk, na.rm = TRUE)) %>%
  ungroup() 

ai_risk_summary_overall <- ai_risk_by_chemgroup %>%
  group_by(year) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE),
            bumblebees_risk = sum(bumblebees_risk, na.rm = TRUE),
            obicornis_risk = sum(obicornis_risk, na.rm = TRUE),
            ocornuta_risk = sum(ocornuta_risk, na.rm = TRUE)) %>%
  ungroup()


total_data <- pivot_longer(ai_risk_summary_overall,
                           cols = !year, names_to = "Taxon group", values_to = "Risk")
total_data$`Taxon group` <- as.character(lapply(total_data$`Taxon group`, function(x) strsplit(x,"_")[[1]][1]))
total_data$`Taxon group` <- factor(total_data$`Taxon group`)
levels(total_data$`Taxon group`) <- c("Bumblebees", "Earthworms", "Honeybees","Lacewings", "Solitary bee O.b.", "Solitary bee O.c.", "Springtails", "Wasps")


colpal <- carto_pal(8, "Safe")

#raw scale
risk_temp <- ggplot(data = total_data) +
  geom_line(aes(x = year, y = Risk, color = `Taxon group`, linetype = `Taxon group`), 
            linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk - absolute")+xlab("Year")+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

ggsave(plot = risk_temp, width = 15, height = 10, units = "cm", dpi = 300,
       filename = "/data/outputs/overall_risk_absolute_temporal.png")

# plot risk on a relative scale (relative to 1994)

ai_risk_summary_relative <- ai_risk_summary_overall
ai_risk_summary_relative[2:9] <- apply(ai_risk_summary_relative[2:9], 2, function(x) x/x[1])

total_data_relative <- pivot_longer(ai_risk_summary_relative,
                                    cols = !year, names_to = "Taxon group", values_to = "Risk")
total_data_relative$`Taxon group` <- as.character(lapply(total_data_relative$`Taxon group`, function(x) strsplit(x,"_")[[1]][1]))
total_data_relative$`Taxon group` <- factor(total_data_relative$`Taxon group`)
levels(total_data_relative$`Taxon group`) <- c("Bumblebees", "Earthworms", "Honeybees","Lacewings", "Solitary bee O.b.", "Solitary bee O.c.", "Springtails", "Wasps")

risk_temp_rel <- ggplot(data = total_data_relative) +
  geom_line(aes(x = year, y = Risk, color = `Taxon group`, linetype = `Taxon group`), 
            linewidth = 1.3) +
  geom_hline(aes(yintercept = 1), linetype = "dotted", color = "black")+
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk - relative to 1994")+xlab("Year")+
  theme_classic()

ggsave(plot = risk_temp_rel, width = 15, height = 10, units = "cm", dpi = 300,
       filename = "/data/outputs/overall_risk_relative_temporal.png")



## How has the contribution of different chemical groups to risk changed over time?

total_data_chemgroup <- pivot_longer(ai_risk_summary,
                                     cols = !c(year, chemgroup), names_to = "Taxon group", values_to = "Risk")
total_data_chemgroup$`Taxon group` <- as.character(lapply(total_data_chemgroup$`Taxon group`, function(x) strsplit(x,"_")[[1]][1]))
total_data_chemgroup$`Taxon group` <- factor(total_data_chemgroup$`Taxon group`)
levels(total_data_chemgroup$`Taxon group`) <- c("Bumblebees", "Earthworms", "Honeybees","Lacewings", "Solitary bee O.b.", "Solitary bee O.c.", "Springtails", "Wasps")



p1 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Honeybees",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to honeybees")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p2 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Bumblebees",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to bumblebees")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p3 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Solitary bee O.b.",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to Osmia bicolor")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p4 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Solitary bee O.c.",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to Osmia cornutus")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()


p5 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Wasps",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to parasitic wasps")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p6 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Lacewings",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to lacewings")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()


p7 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Earthworms",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to earthworms")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  theme_classic()

p8 <- ggplot(data = total_data_chemgroup[total_data_chemgroup$`Taxon group` == "Springtails",]) +
  geom_line(aes(x = year, y = Risk, color = chemgroup,
                linetype = chemgroup), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to springtails")+xlab("Year")+
  guides(color=guide_legend(title="Target group"),
         linetype = guide_legend(title = "Target group"))+
  theme_classic()

risk_chemgroup <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8, labels="AUTO",
          nrow = 4,ncol= 2, common.legend = TRUE)


ggsave(plot = risk_chemgroup, width = 16, height = 29, units = "cm", dpi = 300,
       filename = "/data/outputs/risk_chemgroup.png")



## risk by type

ai_risk_summary_type <- ai_risk_by_type %>%
  group_by(year, type_slct) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE),
            bumblebees_risk = sum(bumblebees_risk, na.rm = TRUE),
            obicornis_risk = sum(obicornis_risk, na.rm = TRUE),
            ocornuta_risk = sum(ocornuta_risk, na.rm = TRUE)) %>%
  ungroup()

total_data_type <- pivot_longer(ai_risk_summary_type,
                                cols = !c(year, type_slct), names_to = "Taxon group", values_to = "Risk")
total_data_type$`Taxon group` <- as.character(lapply(total_data_type$`Taxon group`, function(x) strsplit(x,"_")[[1]][1]))
total_data_type$`Taxon group` <- factor(total_data_type$`Taxon group`)
levels(total_data_type$`Taxon group`) <- c("Bumblebees", "Earthworms", "Honeybees","Lacewings", "Solitary bee O.b.", "Solitary bee O.c.", "Springtails", "Wasps")
total_data_type$`Active ingredient class` <- as.character(lapply(total_data_type$type_slct, function(x) strsplit(x,"_")[[1]][1]))
total_data_type$`Target group` <- as.character(lapply(total_data_type$type_slct, function(x) strsplit(x,"_")[[1]][2]))


colpal <- c(carto_pal(13, "Safe"),carto_pal(12, "Safe"))


p1.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Honeybees" & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to honeybees\nby insecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()


p2.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Bumblebees" & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to bumblebees\nby insecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p3.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Solitary bee O.b." & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to Osmia bicornis\nby insecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p4.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Solitary bee O.c." & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to Osmia cornutus\nby insecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()



p5.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Lacewings" & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to lacewings\nby insecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()


p6.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Wasps" & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = c(colpal, "black")) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to wasps by\ninsecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

p7.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Earthworms" & total_data_type$`Target group` == "Fungicide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = colpal) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to earthworms\nby fungicide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  theme_classic()


p8.1 <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Springtails" & total_data_type$`Target group` == "Insecticide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = c(colpal)) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to springtails\nby insecticide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  theme_classic()



risk_chemtype <- ggarrange(p1.1,p2.1,p3.1,p4.1,p5.1,p6.1,p7.1,p8.1, 
                           labels="AUTO", nrow = 4,ncol= 2)


ggsave(plot = risk_chemtype, width = 23, height = 29, units = "cm", dpi = 300,
       filename = "/data/outputs/risk_chemtype.png")


wasps_fungi <- ggplot(data = total_data_type[total_data_type$`Taxon group` == "Wasps" & total_data_type$`Target group` == "Fungicide",]) +
  geom_line(aes(x = year, y = Risk, color = `Active ingredient class`,
                linetype = `Active ingredient class`), linewidth = 1.3) +
  scale_color_manual(values = c(colpal, "black")) +
  #scale_y_continuous(trans = "log1p",
  #                   labels = label_number(suffix = " K", scale = 1e-3)) +
  ylab("Risk to wasps by\nfungicide class")+xlab("Year")+
  guides(color=guide_legend(title="Active ingredient\nclass"),
         linetype = guide_legend(title="Active ingredient\nclass"))+
  scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6)) +
  theme_classic()

ggsave(plot = wasps_fungi, width = 10, height = 8, units = "cm", dpi = 300,
       filename = "/data/outputs/wasps_fungicides.png")


## How has the spatial pattern of risk changed over time?

ai_1994 <- ai_risk_overall[ai_risk_overall$year == 1994,]
ai_2016 <- ai_risk_overall[ai_risk_overall$year == 2016,]

ai_change <- cbind(ai_2016[,2:4],ai_2016[,5:12]-ai_1994[,5:12])

p1.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = honeybees_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right") +
  ggtitle("Honeybees")

p2.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = bumblebees_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Bumblebees")


p3.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = obicornis_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Osmia bicornis")

p4.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = ocornuta_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Osmia cornuta")


p5.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = lacewings_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  #xlim(100000, 700000) +
  #ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Lacewings")

p6.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = wasps_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Parasitic wasps")


p7.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = earthworms_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Earthworms")

p8.2 <- ggplot() +
  #geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_change, 
            aes(x = E, y = N, fill = springtails_risk)) +
  scale_fill_continuous(type = "viridis",
                        name = "") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  theme_map() +
  theme(legend.position = "right")+
  ggtitle("Springtails")

risk_change <- ggarrange(p1.2,p2.2,p3.2,p4.2,p5.2,p6.2,p7.2,p8.2, 
                           labels="AUTO", nrow = 4,ncol= 2) 

risk_change <- annotate_figure(
  risk_change,
  top = text_grob("Change in exposure risk 1994-2016",
                  face = "bold", size = 14))

ggsave(plot = risk_change, width = 23, height = 29, units = "cm", dpi = 300,
       filename = "/data/outputs/risk_change.png")
