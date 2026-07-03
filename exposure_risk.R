#######################################################
## Risk calculations from pesticide exposure
## Author: Francesca Mancini
## Date created: 2023-09-06
#######################################################

# devtools::install_github("colinharrower/BRCmap")
library(tidyverse)
library(BRCmap)
library(ggplot2)
library(ggthemes)
library(terra)
library(tidyterra)
library(ggpubr)
library(scales)
library(labelr)

inputs_path <- "/data/inputs/"
outputs_path <- "/data/outputs/"

# read weight applied data
ai_weight <- readr::read_csv(file.path(inputs_path, "Weight_applied.csv"))

# read toxicity data
tox_data <- read.csv(file.path(inputs_path, "toxicity_data_updated2_MG update.csv"),
                     na.strings = "*") %>%
  mutate(ai = tolower(ai)) 


regions <- unique(ai_weight$region)
pesticide_list <- unique(ai_weight$ai)

length(pesticide_list)
# there are 276 active ingredients in the data
# we need to filter to only those present in the toxicity data

ai_weight <- ai_weight %>%
  mutate(ai = tolower(ai)) # first uncapitalise all


ai_weight_filt <- ai_weight %>% # then filter
  filter(ai %in% tox_data$ai)


length(unique(ai_weight_filt$ai))
# only 179, we are missing some active ingredients in the data


all_ai <- unique(ai_weight$ai)

tox_data %>%
  filter(!(ai %in% all_ai))

# ai    application chemgroup
# 1 2-methoxyethylmercury acetate Seed treatment Fungicide
# 2           mefentrifluconazole          Spray Fungicide
# 3                   chlormequat          Spray Herbicide
# 4                      ethephon          Spray Herbicide
# 5              trinexapac-ethyl          Spray Herbicide
# 6                   cinmethylin          Spray Herbicide

# the above 6 chemicals with toxicity data are not in the PUS data

rm(ai_weight)

# summarise weight applied by site and year and plot
ai_weight_summary <- ai_weight_filt %>%
  group_by(year, gr) %>%
  summarise(weight_applied = sum(weight_applied)) %>%
  left_join(unique(select(ai_weight_filt, gr, E, N)))

ggplot() +
  geom_path(data = UK, aes(x = long, y = lat, group = group)) +
  geom_tile(data = ai_weight_summary, 
            aes(x = E, y = N, fill = weight_applied)) +
  scale_fill_continuous(type = "viridis", 
                        name = "Weight applied") +
  # xlim(100000, 700000) +
  # ylim(0, 700000) +
  coord_equal() +
  facet_wrap(~year) +
  # ggtitle("Exposure to Zeta-cypermethrin") +
  # theme_map() +
  theme(legend.position = "bottom")

# summarise across all sites and years
# we used these summaries to pick the main substance types to produce outputs for
ai_weight_summary_all_years <- ai_weight_filt %>%
  group_by(ai) %>%
  summarise(weight_applied = sum(weight_applied)) %>%
  left_join(select(tox_data, ai, chemgroup, type))

write.csv(ai_weight_summary_all_years, 
          "./ai_weight_summary.csv",
          row.names = FALSE)

type_weight_summary_all_years <- ai_weight_filt %>%
  group_by(ai) %>%
  summarise(weight_applied = sum(weight_applied)) %>%
  left_join(select(tox_data, ai, chemgroup, type)) %>%
  group_by(type) %>%
  summarise(weight_applied = sum(weight_applied))

write.csv(type_weight_summary_all_years, 
          "./type_weight_summary.csv",
          row.names = FALSE)


# join the pesticide data with the toxicity data
ai_weight_tox <- left_join(select(ai_weight_filt, -chemgroup), 
                           tox_data, by = "ai") 

# calculate the risk metrics

ai_risk <- ai_weight_tox %>%
  mutate(weight_pha_5km = weight_applied/2500) %>% # first calculate the application rate in kg/ha
  mutate(earthworms_risk = (weight_pha_5km * 1.33)/earthworms_lwst_PNEC,
         springtails_risk = (weight_pha_5km * 1.33)/fcandida_lwst_PNEC,
         lacewings_risk = (weight_pha_5km * 1000)/lacewings_lwst,
         wasps_risk = (weight_pha_5km * 1000)/wasps_lwst,
         hb_oral_dose = case_when(
           application == "Spray" ~ weight_pha_5km * 6.4,
           application == "Seed treatment" ~ weight_pha_5km * 1.08,
           application == "Soil application" ~ weight_pha_5km * 2.99),
         hb_contact_dose = weight_pha_5km * (0.0114*1000),
         bb_oral_dose = case_when(
           application == "Spray" ~ weight_pha_5km * 10,
           application == "Seed treatment" ~ weight_pha_5km * 1.67,
           application == "Soil application" ~ weight_pha_5km * 2.99),
         bb_contact_dose = weight_pha_5km * (0.0146*1000),
         sb_oral_dose = case_when(
           application == "Spray" ~ weight_pha_5km * 0.7,
           application == "Seed treatment" ~ weight_pha * 0.12,
           application == "Soil application" ~ weight_pha_5km * 2.99),
         sb_contact_dose = weight_pha_5km * (0.00184*1000)) %>%
  mutate(honeybees_exp = (hb_oral_dose + hb_contact_dose)*2,
         bumblebees_exp = (bb_oral_dose + bb_contact_dose)*2,
         solbees_exp = (sb_oral_dose + sb_contact_dose)*2) %>%
  mutate(honeybees_risk = honeybees_exp/honeybees_lwst_PNEC,
         bumblebees_risk = bumblebees_exp/(honeybees_lwst_PNEC/2.75),
         obicornis_risk = solbees_exp/(honeybees_lwst_PNEC/0.84),
         ocornuta_risk = solbees_exp/(honeybees_lwst_PNEC/1.795)) %>%
  select(year, E, N, gr, crop_type, crop_area, chemgroup, type, ai,
         treatment, weight_pha_5km, earthworms_risk,
         springtails_risk, lacewings_risk, wasps_risk,
         honeybees_risk, bumblebees_risk, obicornis_risk, ocornuta_risk)




saveRDS(ai_risk, "./ai_risk_corrected.rds")

# for how many active ingredients could we calculate risk metrics?

ai_risk_summary <- ai_risk %>%
  group_by(ai) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE)) %>%
  ungroup() 

ai_by_group <- ai_risk_summary %>%
  summarise(earthworms_risk = sum(earthworms_risk > 0),
            springtails_risk = sum(springtails_risk > 0),
            lacewings_risk = sum(lacewings_risk > 0),
            wasps_risk = sum(wasps_risk > 0),
            honeybees_risk = sum(honeybees_risk > 0))


ai_final_list <- select(ai_risk_summary, ai)

write.csv(ai_final_list, "./ai_final_list.csv", row.names = FALSE)
write.csv(ai_by_group, "./ai_by_group.csv", row.names = FALSE)

ai_lacewings <- ai_risk_summary %>%
  select(ai, lacewings_risk) %>%
  filter(lacewings_risk >0) %>%
  select(ai)

tox_data_lacewings <- tox_data %>%
  select(ai, lacewings_lwst) %>%
  filter(!(is.na(lacewings_lwst))) %>%
  select(ai)


tox_data_lacewings %>%
  filter(!(ai %in% ai_lacewings$ai))



ai_springtails <- ai_risk_summary %>%
  select(ai, springtails_risk) %>%
  filter(springtails_risk >0) %>%
  select(ai)

tox_data_springtails <- tox_data %>%
  select(ai, fcandida_lwst_PNEC) %>%
  filter(!(is.na(fcandida_lwst_PNEC))) %>%
  select(ai)

  
tox_data_springtails %>%
  filter(!(ai %in% ai_springtails$ai))


# summarise by main chemical group 
# (insecticide, fungicide, herbicide and molluscicide)

ai_risk_by_chemgroup <- ai_risk %>%
  group_by(year, gr, E, N, chemgroup) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE),
            bumblebees_risk = sum(bumblebees_risk, na.rm = TRUE),
            obicornis_risk = sum(obicornis_risk, na.rm = TRUE),
            ocornuta_risk = sum(ocornuta_risk, na.rm = TRUE)) %>%
  ungroup() 

# summarise across all active ingredients

ai_risk_overall <- ai_risk %>%
  group_by(year, gr, E, N) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE),
            bumblebees_risk = sum(bumblebees_risk, na.rm = TRUE),
            obicornis_risk = sum(obicornis_risk, na.rm = TRUE),
            ocornuta_risk = sum(ocornuta_risk, na.rm = TRUE)) %>%
  ungroup() 


saveRDS(ai_risk_by_chemgroup, "./ai_risk_chemgroup_corrected.rds")
saveRDS(ai_risk_overall,"./ai_risk_overall_corrected.rds")


# summarise across substance type
# as there are more than 70 substance types we picked 
# some major groups based on total weight applied across all years and crops 

chemtypes_slct <- c(
  "Chloronitrile",
  "Triazole" ,
  "Carbamate",
  "Morpholine",
  "Strobilurin", 
  "Triazolinthione",
  "Benzimidazole",
  "Urea",
  "Dinitroaniline",
  "Organophosphate",
  "Thiocarbamate",
  "Aryloxyalkanoic acid",
  "Chloroacetamide",
  "Benzamide",
  "Triazine",
  "Oxyacetamide",
  "Pyridine compound",
  "Triazinone",
  "Pyrethroid",
  "Neonicotinoid",
  "Organochlorine")

ai_risk_by_type <- ai_risk %>%
  mutate(type_slct = case_when(
    type %in% chemtypes_slct ~ paste(type, chemgroup, sep = "_"),
    TRUE ~ paste("other", chemgroup, sep = "_"))) %>%
  mutate(type_slct = case_when(
    type_slct == "Benzamide_Fungicide" ~ "other_Fungicide", # there is a benxamide fungicide, which has not been selected, so we want to change this to "other"
    TRUE ~ type_slct)) %>% 
  group_by(year, gr, E, N, type_slct) %>%
  summarise(earthworms_risk = sum(earthworms_risk, na.rm = TRUE),
            springtails_risk = sum(springtails_risk, na.rm = TRUE),
            lacewings_risk = sum(lacewings_risk, na.rm = TRUE),
            wasps_risk = sum(wasps_risk, na.rm = TRUE),
            honeybees_risk = sum(honeybees_risk, na.rm = TRUE),
            bumblebees_risk = sum(bumblebees_risk, na.rm = TRUE),
            obicornis_risk = sum(obicornis_risk, na.rm = TRUE),
            ocornuta_risk = sum(ocornuta_risk, na.rm = TRUE)) %>%
  ungroup() 

saveRDS(ai_risk_by_type, "./ai_risk_type_corrected.rds")

## create and save raster stacks ----

create_raster_stack <- function(df, years, column, chemgroup, out_path, overwrite){

  ras_stack <- rast(
    lapply(years, function(x) 
      rast(df %>%
             filter(year == x) %>%
             select(E, N, column) ,
           type = "xyz",
           crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +nadgrids=OSTN15_NTv2_OSGBtoETRS.gsb +units=m +no_defs +type=crs",
           extent = ext(c(127500, 657500, 7500, 657500))
      )
    ))
  
  
  names(ras_stack) <- paste0("y_", years)
  
  writeRaster(ras_stack, 
              file.path(out_path, 
                        paste0(column, "_", 
                               ifelse(!is.null(chemgroup), 
                                      chemgroup,
                                      "allchem"), ".tif")),
              overwrite = overwrite)
  
}

# create layers of risk from all chemicals for each invertebrate group
columns <- names(ai_risk_overall)[5:12]

for(c in 1:length(columns)){
  
  df <- ai_risk_overall %>%
    select(E, N, columns[c], year)
  
  create_raster_stack(df = df,
                      years = unique(df$year),
                      column = columns[c],
                      chemgroup = NULL, 
                      out_path = "/data/outputs",
                      overwrite = TRUE)
}

# create layers of risk from major chemical groups for each invertebrate group

chemgroups <- unique(ai_risk_by_chemgroup$chemgroup)

for(i in 1:length(chemgroups)){
  
  for(j in 1:length(columns)){
    
    df <- ai_risk_by_chemgroup %>%
      filter(chemgroup == chemgroups[i]) %>%
      select(E, N, columns[j], year)
    
  create_raster_stack(df = df,
                      years = unique(df$year),
                      column = columns[j],
                      chemgroup = chemgroups[i], 
                      out_path = "/data/outputs",
                      overwrite = TRUE)
  }
}

# create layers of risk from major substance types for each invertebrate group

chemtypes <- unique(ai_risk_by_type$type_slct)

single_ai <- c( # single active ingredient, exclude for licensing issues
  "Chloronitrile",
  "Organochlorine",
  "Oxyacetamide",
  "Triazolinthione"
)

exclude <- gremlr(single_ai, chemtypes, vals = TRUE)


chemtypes <- chemtypes[!(chemtypes %in% exclude)]


for(i in 1:length(chemtypes)){
  
  for(j in 1:length(columns)){
    
    df <- ai_risk_by_type %>%
      filter(type_slct == chemtypes[i]) %>%
      select(E, N, columns[j], year)
    
    create_raster_stack(df = df,
                        years = unique(df$year),
                        column = columns[j],
                        chemgroup = chemtypes[i], 
                        out_path = "/data/outputs",
                        overwrite = TRUE)
  }
}


