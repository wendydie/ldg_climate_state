# 0.1 data downloading from PBDB
url <- modify_url(base_url, query = query)
response <- GET(url)
content_raw <- content(response, "text")
data <- read.csv(text = content_raw, skip = 21, header=TRUE)

# 0.2 data cleaning                      "lat", "lng", "max_ma", "min_ma", "early_interval", "late_interval",
data_clean <- data %>%
  select(all_of(selected_columns))%>%
  filter(!is.na(occurrence_no), !is.na(reference_no), !is.na(collection_no),
         !is.na(genus), 
         !is.na(lat), !is.na(lng), 
         !is.na(max_ma),!is.na(min_ma),
         lat != "",lng != "",genus != "",max_ma != "",min_ma != ""
         )

# Add a column to flag uncertain records
data_clean <- data_clean %>%
  mutate(is_uncertain = str_detect(genus, "\\s(sp\\.|cf\\.|aff\\.|indet|nov\\.sp\\.|s\\.l\\.|s\\.s\\.|morphotype.*)"))

# Filter data based on certainty (optional)
data_clean <- data_clean %>%
  filter(!is_uncertain) # Remove uncertain records if not needed

# Remove suffixes from genus names
data_clean <- data_clean %>%
  mutate(genus = str_remove(genus, "\\s(indet\\.|sp\\.|cf\\.|aff\\.|nov\\.sp\\.|sp\\.nov\\.|s\\.l\\.|s\\.s\\.|spp\\.|ex\\sgr\\.|incertae\\ssedis|morphotype.*)$"))%>%
  mutate(genus = str_extract(genus, "^[^\\(]+")) %>% # Extract genus name before the parenthesis
  mutate(genus = str_trim(genus)) # Trim leading/trailing spaces

# Calculate genus frequencies and split by the first letter
genus_freq <- data_clean %>%
  count(genus, sort = TRUE) %>%
  # filter(n <= 2) %>% # Focus only on low-frequency genera
  mutate(first_letter = substr(genus, 1, 1)) %>%
  arrange(first_letter, genus)

# Correct genus names by comparing with neighbors
genus_freq <- genus_freq %>%
  group_by(first_letter) %>%
  mutate(
    corrected_genus = sapply(
      seq_along(genus),
      function(i) {
        current <- genus[i]
        neighbors <- genus[max(1, i - 1):min(n(), i + 1)] # Get neighbors
        distances <- stringdist(current, neighbors)
        closest <- neighbors[which.min(distances)]
        if (min(distances) == 1) closest else current
      }
    )
  )%>%
  ungroup()

# Merge corrected genus names back to the original data
data_clean <- data_clean %>%
  left_join(genus_freq %>% select(genus, corrected_genus), by = "genus") %>%
  mutate(corrected_genus = ifelse(is.na(corrected_genus), genus, corrected_genus))

# 0.3 To-Do: Paleolatitude Reconstruction Annotation
#######################
######################
#######################
######################
#######################

# 0.4 Prepare data for calculating Richness in stages
data_clean <- data_clean %>%
  filter(paleolat2 <= 90 & paleolat2 >= -90)
data_clean <- data_clean %>%
  filter(paleolng2 <= 180 & paleolng2 >= -180)
data_clean['mean_ma'] <- (data_clean['max_ma'] + data_clean['min_ma']) / 2
stages <- read.csv("./data/stages.csv")
stage_stbin <- stages[,c('stage', 'stbin', 'short', 'bottom', 'mid', 'top',
                         'dur')]
data_clean <- data_clean %>%
  rowwise() %>%
  mutate(
    stage = ifelse(length(stage_stbin$stage[which(mean_ma <= stage_stbin$bottom & mean_ma > stage_stbin$top)]) > 0,
                   stage_stbin$stage[which(mean_ma <= stage_stbin$bottom & mean_ma > stage_stbin$top)][1], NA),
    stbin = ifelse(length(stage_stbin$stbin[which(mean_ma <= stage_stbin$bottom & mean_ma > stage_stbin$top)]) > 0,
                   stage_stbin$stbin[which(mean_ma <= stage_stbin$bottom & mean_ma > stage_stbin$top)][1], NA),
    stage_mid = ifelse(length(stage_stbin$mid[which(mean_ma <= stage_stbin$bottom & mean_ma > stage_stbin$top)]) > 0,
                       stage_stbin$mid[which(mean_ma <= stage_stbin$bottom & mean_ma > stage_stbin$top)][1], NA)
  ) %>%
  ungroup()
saveRDS(data_clean, file = "data/data_clean.rds")
# saveRDS(data, file = "data/data.rds")
# data_clean <- readRDS("data/data_clean.rds")
