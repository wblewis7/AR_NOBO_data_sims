# AR_NOBO_data_sims
Data and code for simulations and empirical analyses of northern bobwhite point count and ARU data from Arkansas, USA. From Lewis, W. B., C. Johnson, and J. A. Martin. Assessing the efficacy and cost-effectiveness of integrating autonomous recording units and point-count surveys for population monitoring of northern bobwhite (Colinus virginianus).
Contact Information: William Lewis, wblewis7@gmail.com, University of Georgia
---

---

# Metadata

# AR_NOBO_data.gzip
Empirical data comes for the northern bobwhite (Colinus virginianus) project are stored in the "AR_NOBO_data" gzip file. Calling bobwhite groups (coveys) were surveyed via point counts and autonomous recording units (ARUs) in the fall of 2022 at eight study sites across
Arkansas, USA. Trained observers performed 600-m radius point counts at 3 - 6 survey locations at each study study (29 total locations), recording the number of calling coveys and the distance to the observer. Distances were grouped into 6 100-m bins for analysis. Surveys were performed from 45 minutes before to 15 minutes after sunrise and were repeated an average of 3 times (range 1 - 5). SM mini Wildlife Acoustic ARUs were deployed at 2 - 5 locations surveyed by point counts in each study site (26 total locations), set to record each day during the same time period as point counts. ARU recordings from days with inclement weather (measured precipitation > 0 or wind speeds > 16 kmph) were removed. Deployment dates for ARUs ranged from October 1 - Novemeber 10. Bobwhite calls were detected in ARU recordings using a previously developed convolutional neural network and a score threshold of 0.95. In a preliminary analysis, the detection radius for ARUs was determined to be smaller than for point counts. Point counts were performed at the same locations as ARUs, so point counts surveyed the entire area surveyed by ARUs plus an additional area.
## y.PC
Formatted point count survey data giving the total number of bobwhite covey detections per survey. Rows represent unique survey locations while columns represent repeat surveys at the same location. NA values represent no survey for that number of visits at a location.
## ydb.PC
Formatted point count survey data giving the number of bobwhite covey detections in each detection bin per survey. The first dimension of the array corresponds to unique survey locations, the second dimension represents the 6 distance bins (0-100m, 100-200m, 200-300m,
300-400m, 400-500m, 500-600m), and the third dimension represents repeat surveys at the same location.
## nareaeffects
The number of study area effects to estimate in the model. This number is the number of study areas minus 1, since the first site is treated as the intercept.
## nB
The number of distance bands from the point count distance sampling.
## db
Break points of distance bands from the point count distance sampling.
## pix
Relative area of each distance band from the point count distance sampling.
## point.area.PC
Area (square meters) surveyed by each point count.
## point.area.ARU
Area (square meters) surveyed by ARU.
## point.area.DIFF
Additional area (square meters) surveyed by point counts but not by ARUs.
## npoints
Number of unique survey locations
## area
Matrix giving the study area for each survey location. Rows represent unique survey locations. Columns represent study area effects estimated in the model. Values of 1 represent which study area effect to use for each survey location. Survey locations at the reference
study area are denoted with all 0s.
## nvisits
The number of repeat point count visits performed at each survey location.
## nareas
The total number of survey areas (8)
## nvalidate
The number of ARU recordings manually validated for false positives.
## FP
The number of detections identified by the automated classifier per validated recording determined to be false positives through manual validation.
## n.ARU
The number of ARUs used in the study
## n.A.times
The number of recordings for each ARU which detected at least one covey call.
## A.times
For each ARU, gives the day of the season for recordings at which at least one covey call was detected. Rows represent ARUs, and columns represent the days of recordings at that ARU which detected calls. The number of non-NA values for each row are indexed by n.A.times.
## v.ARU
Formatted ARU data giving the total number of bobwhite covey calls detected per recording. Rows represent ARUs and columns represent days of recordings during the deployment period. NA values represent days without recordings.
## y.ARU
Formatted ARU data giving the presence or absence of bobwhite covey calls detected on recordings. Rows represent ARUs and columns represent days of recordings during the deployment period. NA values represent days without recordings.
## ARU.pointID
Vector giving the point location at which each ARU was deployed.
## ARU.noise
Matrix giving standardized values of background noise on ARU recordings. Rows represent ARUs and columns represent days of recordings during the deployment period. Background noise was used to assess variation in the calling/detection parameter (delta)
## maxdates
The maximum number of days of recording across ARUs.
 

<br />
<br />

# Analyzing_AR_NOBO_data.R
Sample code for jointly estimating abundance of bobwhite coveys from point count and ARU data is provided in Analyzing_AR_NOBO_data.R. Abundance was modeled as varying by study area, with seperate observation processes for point counts and ARUs. The statistical model
follows that of Doser et al. (2021). Integrating automated acoustic vocalization data and point count surveys for estimation of bird abundance, with a few exceptions. Point counts survey a larger area than ARUs based on a preliminary analysis; counts were performed right at ARU locations, so point counts surveyed the entire ARU area plus an additional area. We modeled abundance separatley for the area surveyed by ARUs and the expanded area between the ARU and point count detection radius using the same parameters but an offset for area. These two abundance estimates were summed to determine abundance within the point count sampling range. We also modified the false positive estimation process, using a Poisson process and the number of false positives per validated recording rather than a
Hypergeometric distribution. Finally, we model variation in the parameter for the average number of calls detected/covey (delta) based on study area and background noise on recordings.
 

<br />
<br />

# Generating_analyzing_NOBO_simulation_data.R
Sample code for simulations to test the efficacy of the Doser et al (2021) model to jointly estimate abundance of bobwhite coveys from point counts and ARUs. Code is provided to simulate point count, ARU, and manually-validated false positive data and then to run the model of Doser et al. (2021) on the simulated data. Data were simulated for a hypothetical 2016 ha study area with average bird density of 0.1, 0.5, 0.9, 1.3, 1.7, 2.1, and 2.5 birds/ha (and assuming 12 birds per covey). Simulated data was generated with surveying either 2, 4, 6, or 8 locations (11 - 56% of study area), either all with point counts and ARUs, all with point counts and half with ARUs, all with point counts and none with ARUs, all with ARUs and half with point counts, all with ARUs and none with point counts, and half with point counts and the other half with ARUs. Areas surveyed by point counts were surveyed either 2, 3, 4, or 5 repeat times and ARUs were deployed for 7, 14, 21, or 28 days.
Simulation parameters were based on estimated values from analysis of the Arkansas empirical data. Analysis of the Arkansas data revealed that many of the parameters for the ARU detection process were not modeled well. We therefore used a slighlty different method of generating ARU data compared to how the data are analyzed in the statistical model. In particular, we increased the variability in the average number of calls detected/covey (delta) to make the data more biologically realistic.
