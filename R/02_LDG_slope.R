# 02 calculating the slope of LDG in time bins (N,S,combined)

# Step 1: Split Assemblage into paleolat2 and paleolng2
combined_rich <- combined_rich %>%
  separate(Assemblage, into = c("paleolat2", "paleolng2"), sep = "_", convert = TRUE)

# Step 2: Calculate combined slope by stbin
combined_slope <- combined_rich %>%
  group_by(stbin) %>%
  summarise(
    combined_slope = coef(lm(qD ~ paleolat2))[2], # Calculate slope for qD ~ paleolat2
    .groups = 'drop'
  )

# Step 3: Calculate southern and northern slopes
southern_slope <- combined_rich %>%
  filter(paleolat2 < 0) %>%
  group_by(stbin) %>%
  summarise(
    southern_slope = coef(lm(qD ~ paleolat2))[2], # Slope for southern hemisphere
    .groups = 'drop'
  )

northern_slope <- combined_rich %>%
  filter(paleolat2 >= 0) %>%
  group_by(stbin) %>%
  summarise(
    northern_slope = coef(lm(qD ~ paleolat2))[2], # Slope for northern hemisphere
    .groups = 'drop'
  )

# Step 4: Combine the results
final_results <- combined_slope %>%
  left_join(southern_slope, by = "stbin") %>%
  left_join(northern_slope, by = "stbin")

# View the final results
print(final_results)
