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
## nB
The number of distance bins from the point count distance sampling.
## db
Break points of distance bins from the point count distance sampling.
## pix
Relative area of each distance bin from the point count distance sampling.
## areaPC
Area (square meters) surveyed by each point count.
## areaARU
Area (square meters) surveyed by ARUs.
## areaDiff
Additional area (square meters) surveyed by point counts but not by ARUs.
## npoints
Number of unique survey locations
## site
Vector giving the study site for each survey location.
## nsites
The number of study sites (8).
## nvisits
The number of repeat point count visits performed at each survey location.
## n.ARU
The number of ARUs used in the study (26).
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
Maximum number of days of recording across ARUs.
## n.ARU.val
Number of recordings from which a subset of ARU detections were manually validated.
## which.v.ARU
Indexing vector giving the ARUs in y.ARU with at least one recording with a detected bobwhite call.
## n.which.v.ARU
Number of ARUs with at least one recording with a detected bobwhite call.
## ARUID.val
Indexing vector giving the ARU in v.ARU associated with each recording with manual validation.
## ARUID.sub.val and ARU.Atimes.ID.val
Indexing vectors giving the delta estimates associated with each recording with manual validation.
## ARUday.val
Indexing vector giving the day of recording in v.ARU associated with each recording with manual validation.
## k.val
Number of true positives within manually validated detections for recordings with validation.
## n.val
Total number of manually validated detections for recordings with validation.

<br />
<br />

# Analyzing_AR_NOBO_data.R
Sample code for jointly estimating abundance of bobwhite coveys from fall distance-sampling point count and ARU data is provided in Analyzing_AR_NOBO_data.R. Point count and ARU data are integrated in a Bayesian hierarchical framework. Both datasets jointly estimate the state (abundance) process but have seperate observation processes. We use a modified version of the statistical model of Doser et al. 2021 (Integrating automated acoustic vocalization data and point count surveys for estimation of bird abundance), with modifications by Nolan et al. 2024 (Effects of management practices on Northern Bobwhite Colinus virginianus density in privately owned working forests across the Southeastern United States). Point counts are incorporated via a hierarchical distance-sampling model (Royle et al. 2004 Modeling abundance effects in distance sampling) with a Half-Normal detection function (Addressing temporal variability in bird calling with design and estimation: A Northern Bobwhite example), modified to account for repeat surveys. ARU data are incorporated via a zero-truncated hurdle model and account for false positives and negatives. We model abundance as varying by study site, while the ARU detection process varies by site and background noise of recordings. A subset of detections are manually validated and incorporated with a Hypergeometric formulation to correct for false positives.

<br />
<br />

# Generating_analyzing_NOBO_simulation_data.R
Sample code for generating and analyzing simulated datasets to assess the efficacy of Bayesian hierarchical models to estimate abundance of bobwhite coveys from point counts and ARUs. Code is provided to simulate point count, ARU, and manually-validated false positive data and then to use models similar to in Analyzing_AR_NOBO_data.R for estimating abundance. Data were simulated for a hypothetical 2016 ha study area with average bird density of 0.1, 0.5, 0.9, 1.3, 1.7, 2.1, and 2.5 birds/ha (and assuming 12 birds per covey). Simulated data was generated with surveying either 22, 34, 45, or 56% of the study area (4, 6, 8, or 10 survey locations) with either point counts (model C), ARUs with manual validation of a subset of detections (model AV), or both (AVC). Point counts were repeated 2 - 4 times, while ARUs were deployed for either 14 or 28 days. 100 datasets were simulated for each of the 308 scenarios.
Simulation parameters were based on analysis of point count and ARU data on bobwhite abundance from Arkansas, USA, 2022. Analysis of the Arkansas data suggested that the number of detections across recordings on the same ARU was highly variable, and this varition was poorly explained by estimates of the parameter (delta), which governs the calling rate and detection probability of bobwhite calls on ARUS. We therefore used a slightly different procedure between the data generation process and the statistical model for calculating delta. We believe that this method of data generation creates more biologically realistic data, thus better assessing the performance of the statistical model on ecological datasets; however, we acknowledge that the simulations are unlikely to be able to recover parameters associated with the ARU detection process. We therefore only assessed model performance of density estimates and the associated parameter β_0.
Data are analyzed via Bayesian hierarchical models based on the models of Doser et al, 2021 (Integrating automated acoustic vocalization data and point count surveys for estimation of bird abundance), with modifications by Nolan et al. 2024 (Effects of management practices on Northern Bobwhite Colinus virginianus density in privately owned working forests across the Southeastern United States). Point count data are incorporated via distance-sampling models modified for repeat surveys, while ARU data are incorporated via a zero-truncated hurdle model to account for false positives and false negatives.
