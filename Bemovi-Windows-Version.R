################################################################################
### bemovi phenotyping of paramecium
################################################################################

### load bemovi

library(bemovi); library(ggplot2); library(dplyr); library(cowplot); library(ggcorrplot)

### project directory 

to.data <- 'YOUR PROJECT DIRECTORY/'

setwd(to.data)

### provide two folders
# 1) Video description
# 2) the folder with the raw videos

video.description.folder <- '0 - video description/'

video.description.file <- 'video_description.txt'

raw.video.folder <- '1 - raw/'

# intermediate directories that are automagically created

particle.data.folder <- "2 - particle data/"
trajectory.data.folder <- "3 - trajectory data/"
temp.overlay.folder <- "4a - temp overlays/"
overlay.folder <- "4 - overlays/"
merged.data.folder <- "5 - merged data/"
ijmacs.folder <- "ijmacs/"

### check that video file names are ok

check_video_file_names(to.data, raw.video.folder, video.description.folder, video.description.file)

## path to image J and particle linker and java 

IJ.path <- 'C:/Users/DeLongLab/Documents/ImageJ/ImageJ.exe'

to.particlelinker <- 'C:/Users/DeLongLab/Documents/VideoAnalysisExample/VideoAnalysis/'

java.path <- "C:/Users/DeLongLab/Documents/ImageJ/jre/bin/javaw.exe"

# amount of memory available to ImageJ and Java

memory.alloc <- c(10000)

# video frame rate

fps <- 18.1 # frames per second

# pixels per micrometer

pixel_to_scale <- 10.8 # pixels per micrometer

# difference image segmentation

difference.lag <- 5

# specify threshold

thresholds <- c(30, 255)

check_threshold_values(to.data, raw.video.folder,
                       ijmacs.folder, 20, difference.lag,
                       thresholds, IJ.path, memory.alloc)

### segment and locate particles

locate_and_measure_particles(to.data, raw.video.folder,
                             particle.data.folder, difference.lag,
                             thresholds, min_size = 5, max_size = 1000,
                             IJ.path, memory.alloc)

### link particles

link_particles(to.data, particle.data.folder,
               trajectory.data.folder, linkrange = 10,
               disp = 30, start_vid = 1, memory = memory.alloc)

# merge data 

merge_data(to.data, particle.data.folder, 
           trajectory.data.folder, video.description.folder,
           video.description.file, merged.data.folder)

#### Processing the data

load(paste0(to.data,merged.data.folder,"Master.RData"))

### filtering the data 

trajectory.data.filtered <- filter_data(trajectory.data,50,8,0.5,10)
trajectory.data.filtered$type <- "filtered data"
trajectory.data$type <- "raw data"

morph_mvt <- summarize_trajectories(trajectory.data.filtered, write=T, calculate.median=F, to.data, merged.data.folder)

morph_mvt

morph_mvt %>% group_by(file) %>% summarise(number = n()) %>% print(n = 140)

### save data

write.csv(morph_mvt, file = 'YOUR FILE NAME.csv')
