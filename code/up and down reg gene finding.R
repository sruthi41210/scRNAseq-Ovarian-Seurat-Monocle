library(dplyr)

# Assuming all_markers was already generated
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(10, wt = avg_log2FC)

# Save as Excel-compatible CSV
write.csv(top_markers, "Zhang_Top_Markers.csv", row.names = FALSE)
