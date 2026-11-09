### Bemovi phenotyping of Paramecium --------------------------------------

library(bemovi)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggbiplot)
library(ggcorrplot)

## Paths -------------------------------------------------------------------

proj_dir <- "YOUR_PROJECT_DIRECTORY"
setwd(proj_dir)

to.data                <- file.path(proj_dir, "")
video.description.folder <- "0 - video description/"
video.description.file   <- "video.description.txt"
raw.video.folder         <- "1 - raw video/"

particle.data.folder   <- "2 - particle data/"
trajectory.data.folder <- "3 - trajectory data/"
temp.overlay.folder    <- "4a - temp overlays/"
overlay.folder         <- "4 - overlays/"
merged.data.folder     <- "5 - merged data/"
ijmacs.folder          <- "ijmacs/"

## ImageJ / Java settings --------------------------------------------------

IJ.path         <- "/Applications/Fiji.app/jars" # Path to ImageJ jars directory
to.particlelinker <- "/Users/dinelkathilakarathnegmail.com/Desktop/Bemovi/" # Path to particle linker directory
memory.alloc    <- 50000
fps             <- 30 # Frames per second
pixel_to_scale  <- 76.5 # Pixel-to-scale conversion factor
difference.lag  <- 5 # Difference lag for particle tracking
thresholds      <- c(24, 255)

## 1. Check file names & thresholds ----------------------------------------

check_video_file_names(
  to.data, raw.video.folder,
  video.description.folder, video.description.file
)

check_threshold_values(to.data, raw.video.folder, 
    ijmacs.folder, 0, difference.lag, thresholds, 
    IJ.path, memory.alloc)

## 2. Segment & locate particles ------------------------------------------

locate_and_measure_particles(to.data, raw.video.folder, 
    particle.data.folder, difference.lag, 
    thresholds, min_size = 5, max_size = 1000, 
    IJ.path, memory.alloc)
#Sys.sleep(60)

## 3. Link particles -------------------------------------------------------

link_particles(to.data, particle.data.folder, 
    trajectory.data.folder, linkrange = 5, 
    disp = 20, start_vid = 1, memory = memory.alloc)
#Sys.sleep(60)


## 4. Merge data -----------------------------------------------------------

merge_data(to.data, particle.data.folder, 
    trajectory.data.folder, video.description.folder, 
    video.description.file, merged.data.folder)

load(paste0(to.data, merged.data.folder, "Master.RData"))

## 5. Filter + summarize trajectories --------------------------------------

trajectory.data.filtered <- filter_data(trajectory.data, 200, 9, 0.5, 2.5)
trajectory.data.filtered$type <- "filtered"
trajectory.data$type          <- "raw"

morph_mvt <- summarize_trajectories(
  trajectory.data.filtered,
  write = TRUE,
  calculate.median = FALSE,
  to.data = to.data,
  merged.data.folder = merged.data.folder
)

## 6. PCA on morphology & movement ----------------------------------------

morph_mvt_forpca <- morph_mvt %>%
  select(
    -c(mean_grey, sd_grey, sd_area, sd_perimeter,
       sd_major, sd_minor, sd_ar, duration,
       N_frames, id)
  )

pca_mvt <- prcomp(
  morph_mvt_forpca[, -1],
  center = TRUE, scale. = TRUE
)

all_plot <- ggbiplot(
  pca_mvt,
  groups  = morph_mvt$file,
  ellipse = TRUE
) +
  theme(legend.position = "none")

#save_plot("Phenotype_all_pca.png", all_plot, base_height = 8)

## 7. PCA on file means ----------------------------------------------------

mean_mvt <- morph_mvt_forpca %>%
  group_by(file = as.factor(file)) %>%
  summarise(across(everything(), mean), .groups = "drop")

mean_pca <- prcomp(
  mean_mvt[, -c(1, 2)],
  center = TRUE, scale. = TRUE
)

mean_plot <- ggbiplot(
  mean_pca,
  labels  = unique(morph_mvt$file),
  groups  = unique(morph_mvt$file)
) +
  theme(legend.position = "none")

#save_plot("Phenotype_mean_pca.png", mean_plot, base_height = 6)

## 8. Correlation matrix ---------------------------------------------------

ggcorrplot(cor(mean_mvt[, -c(1, 2)]))

## 9. Example boxplots & scatter ------------------------------------------

ggplot(morph_mvt, aes(x = reorder(file, mean_area), y = mean_area,
                      colour = reorder(file, mean_area))) +
  geom_boxplot()

ggplot(morph_mvt, aes(x = reorder(file, mean_perimeter), y = mean_area,
                      colour = reorder(file, mean_area))) +
  geom_boxplot()

ggplot(mean_mvt, aes(x = mean_ar, y = gross_speed)) +
  geom_point()

## 10. Export data ---------------------------------------------------------

write.csv(morph_mvt, "Paramecium_Phenotype_Data.csv", row.names = FALSE)

