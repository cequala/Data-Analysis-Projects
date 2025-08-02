\[VT\] Find genes related to long survival in GC patients
================
Jongwu Kim
2025-08-02

In baseline timepoint, which means the timepoint before any treatment,
found (ATM+ARID1A+TP53) mutation in SNV/InDel are predictor of long
survival in gastric cancer patients.

``` r
# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# BiocManager::install("limma"); BiocManager::install("preprocessCore"); BiocManager::install("sva")
library(ggplot2); library(FactoMineR); library(ggrepel); library(patchwork); library(reshape2)
library(gridExtra); library(MASS); library(skimr); library(dplyr); library(tidyr); library(openxlsx); library(readr) # For loading data
library(survival); library(survminer); library(maxstat); library(knitr) # For survivial analysis

# Set seed for reproducibility
set.seed(42)
# Use the normal numerical representation
options(scipen = 999)
```

# Prepare Data

``` r
# Load 'VIKTORY-1_best_of_response_for_analysis.xlsx'
br_df <- read.xlsx(file.path(PJT_DIR, "VIKTORY-1_best_of_response_for_analysis.xlsx"), sheet = "VIKTORY-1")
# Check the structure of the data
head(br_df[, 1:10])  # Display the first 10 columns of the data
```

    ##   No. Subject.Enrollment.Number Best.Response.(Cycle)     HER2  EBV    MLH
    ## 1   1                B11-A-1-01                    PD        0  pos Intact
    ## 2   2                B11-A-1-02              SCR Fail        0  neg Intact
    ## 3   3                B11-A-1-03              SCR Fail        0  neg Intact
    ## 4   4                B11-A-1-04                    PD negative <NA> Intact
    ## 5   5                B11-A-1-05                    PD        0  neg Intact
    ## 6   6                B11-A-1-06                    PD       1+  neg Intact
    ##   PD-L1(28-8) PD-L1(22C3) CLDN18  TMB
    ## 1           5          20   <NA>  Low
    ## 2           0           1   <NA>  Low
    ## 3          NA           2   <NA> <NA>
    ## 4           1          NA   <NA>  Low
    ## 5           3           5   <NA>  Low
    ## 6          35           5   <NA> High

``` r
dim(br_df)  # Check the dimensions of the data
```

    ## [1] 34 60

``` r
colnames(br_df)
```

    ##  [1] "No."                                  
    ##  [2] "Subject.Enrollment.Number"            
    ##  [3] "Best.Response.(Cycle)"                
    ##  [4] "HER2"                                 
    ##  [5] "EBV"                                  
    ##  [6] "MLH"                                  
    ##  [7] "PD-L1(28-8)"                          
    ##  [8] "PD-L1(22C3)"                          
    ##  [9] "CLDN18"                               
    ## [10] "TMB"                                  
    ## [11] "MSI"                                  
    ## [12] "MMRD"                                 
    ## [13] "SNV/Indel"                            
    ## [14] "CNV"                                  
    ## [15] "Fusion"                               
    ## [16] "age"                                  
    ## [17] "sex"                                  
    ## [18] "status"                               
    ## [19] "reason.for.end.of.administration"     
    ## [20] "Subject.Enrollment.Number"            
    ## [21] "Subject.Initial"                      
    ## [22] "Base.line.TL.Sum.(cm)"                
    ## [23] "Max.change.from.baseline.(Percentage)"
    ## [24] "Last.Date.of.IP.Administration"       
    ## [25] "Durvalumab.(fix.dose)"                
    ## [26] "First.dose.(Ceralasertib)"            
    ## [27] "Current.Dose.(Ceralasertib)"          
    ## [28] "DR"                                   
    ## [29] "DR.Cycle"                             
    ## [30] "DR.date"                              
    ## [31] "Reason.of.Dose.reduction"             
    ## [32] "progression.1=event.0=none"           
    ## [33] "Start.Date.of.IP.Administration"      
    ## [34] "date.of.progression"                  
    ## [35] "death.or.alive"                       
    ## [36] "date.of.death.or.last.FU"             
    ## [37] "Basline.Biopsy"                       
    ## [38] "Post-Biopsy.(C1D8.pre)"               
    ## [39] "Post-Biopsy.(C2D1.pre)"               
    ## [40] "Basline"                              
    ## [41] "Post-Biopsy.(C1D8.pre)"               
    ## [42] "Post-Biopsy.(C2D1.pre)"               
    ## [43] "BL"                                   
    ## [44] "FU1"                                  
    ## [45] "FU2"                                  
    ## [46] "BL"                                   
    ## [47] "FU1"                                  
    ## [48] "FU2"                                  
    ## [49] "FU3"                                  
    ## [50] "FU4"                                  
    ## [51] "FU5"                                  
    ## [52] "FU6"                                  
    ## [53] "FU7"                                  
    ## [54] "FU8"                                  
    ## [55] "FU9"                                  
    ## [56] "FU10"                                 
    ## [57] "FU11"                                 
    ## [58] "FU12"                                 
    ## [59] "FU13"                                 
    ## [60] "FU14"

``` r
table(br_df$`Best.Response.(Cycle)`)
```

    ## 
    ##         PD    PR (C3)    PR (C5)    PR (C6)     PR(C3)   SCR Fail         SD 
    ##         12          5          2          1          1          4          7 
    ## withdrawal 
    ##          2

``` r
# Drop 'SCR Fail' rows from br_df
br_df <- br_df[br_df$`Best.Response.(Cycle)` != "SCR Fail", ]
nrow(br_df)  # Check the number of rows after dropping 'SCR Fail'
```

    ## [1] 30

``` r
# Drop 'withdrawal' rows from br_df
br_df <- br_df[br_df$`Best.Response.(Cycle)` != "withdrawal", ]
nrow(br_df)  # Check the number of rows after dropping 'SCR Fail'
```

    ## [1] 28

``` r
# Check if there are any NA in `progression.1=event.0=none`
any(is.na(br_df$`progression.1=event.0=none`))  # Check for NA values in the progression column
```

    ## [1] FALSE

``` r
# Rename the `progression.1=event.0=none` column to `PD`
# Ensure column names are unique
colnames(br_df) <- make.unique(colnames(br_df))

# Rename the `progression.1=event.0=none` column to `PD`
br_df <- br_df %>% rename(PD = `progression.1=event.0=none`)
table(br_df$PD)  # Check the distribution of the PD column
```

    ## 
    ##  0  1 
    ## 12 16

``` r
any(is.na(br_df$`death.or.alive`))  # Check for NA values in the death.or.alive column
```

    ## [1] FALSE

``` r
# One-hot encode the "death.or.alive" column
br_df <- br_df |>
  mutate(death.or.alive = ifelse(`death.or.alive` == "Alive", 0, 1)) |>
  rename(Death = `death.or.alive`)
table(br_df$Death)  # Check the distribution of the Death column
```

    ## 
    ##  0  1 
    ## 16 12

``` r
# See the PFS column: "date.of.progression"
br_df$`date.of.progression`[1:5]  # Display the first 5 entries of the PFS column
```

    ## [1] "2024-04-25" "2024-06-12" "2024-07-17" "2024-07-03" "2025-07-23"

``` r
# See the OS column: "date.of.death.or.last.FU"
br_df$`date.of.death.or.last.FU`[1:5]  # Display the first 5 entries of the OS column
```

    ## [1] "2024-08-17" "2024-07-05" "2024-08-15" "2024-07-28" "2025-07-23"

``` r
# See the Basline.Biopsy column: "Basline.Biopsy"
br_df$`Basline.Biopsy`[1:5]  # Display the first 5 entries of the Baseline.Biopsy column
```

    ## [1] 45341 45406 45425 45440 45439

``` r
br_df[c("Subject.Enrollment.Number", "Start.Date.of.IP.Administration", "date.of.progression", "date.of.death.or.last.FU")]
```

    ##    Subject.Enrollment.Number Start.Date.of.IP.Administration
    ## 1                 B11-A-1-01                           45350
    ## 4                 B11-A-1-04                           45419
    ## 5                 B11-A-1-05                           45434
    ## 6                 B11-A-1-06                           45441
    ## 7                 B11-A-1-07                           45447
    ## 9                 B11-A-1-09                           45455
    ## 10                B11-A-1-10                           45455
    ## 11                B11-A-1-11                           45469
    ## 12                B11-A-1-12                           45475
    ## 13                B11-A-1-13                           45468
    ## 14                B11-A-1-14                           45490
    ## 17                B11-A-1-17                           45489
    ## 18                B11-A-1-18                           45497
    ## 19                B11-A-1-19                           45545
    ## 20                B11-A-1-20                           45609
    ## 21                B11-A-1-21                           45623
    ## 22                B11-A-1-22                           45624
    ## 23                B11-A-1-23                           45623
    ## 24                B11-A-1-24                           45637
    ## 26                B11-A-1-26                           45645
    ## 27                B11-A-1-27                           45650
    ## 28                B11-A-1-28                           45672
    ## 29                B11-A-1-29                           45692
    ## 30                B11-A-1-30                           45699
    ## 31                B11-A-1-31                           45700
    ## 32                B11-A-1-32                           45691
    ## 33                B11-A-1-33                           45700
    ## 34                B11-A-1-34                           45741
    ##    date.of.progression date.of.death.or.last.FU
    ## 1           2024-04-25               2024-08-17
    ## 4           2024-06-12               2024-07-05
    ## 5           2024-07-17               2024-08-15
    ## 6           2024-07-03               2024-07-28
    ## 7           2025-07-23               2025-07-23
    ## 9           2025-07-23               2025-07-23
    ## 10          2024-10-01               2025-01-06
    ## 11          2025-01-27               2025-02-17
    ## 12          2025-07-23               2025-07-23
    ## 13          2024-08-06               2024-09-08
    ## 14          2024-09-25               2025-01-09
    ## 17          2024-09-04               2025-01-08
    ## 18          2025-07-23               2025-07-23
    ## 19          2025-07-23               2025-07-23
    ## 20               45713               2025-07-23
    ## 21          2025-07-23               2025-07-23
    ## 22               45680               2025-07-23
    ## 23               45789               2025-07-23
    ## 24               45728               2025-04-12
    ## 26               45785               2025-07-23
    ## 27               45707               2025-05-19
    ## 28               45733               2025-07-23
    ## 29          2025-07-23               2025-07-23
    ## 30          2025-07-23               2025-07-23
    ## 31               45743               2025-04-02
    ## 32          2025-07-23               2025-07-23
    ## 33          2025-07-23               2025-07-23
    ## 34          2025-07-23               2025-07-23

``` r
# Function to convert Excel dates (either as numbers or formatted strings)
convert_excel_date <- function(date_col) {
  # Convert to character first to ensure consistent handling
  date_col <- as.character(date_col)
  
  # Create a result vector
  result <- rep(as.Date(NA), length(date_col))
  
  for (i in seq_along(date_col)) {
    # Skip NA values
    if (is.na(date_col[i])) next
    
    # Try to convert as a regular date string
    date_val <- try(as.Date(date_col[i], format = "%Y-%m-%d"), silent = TRUE)
    
    # If that fails, try as an Excel numeric date
    if (inherits(date_val, "try-error") || is.na(date_val)) {
      # Check if it's numeric
      if (!is.na(suppressWarnings(as.numeric(date_col[i])))) {
        # Convert Excel numeric date (origin is 1899-12-30 for Excel)
        date_val <- as.Date(as.numeric(date_col[i]), origin = "1899-12-30")
      }
    }
    
    result[i] <- date_val
  }
  
  return(result)
}

# Convert the date-related columns to Date format
br_df$`Start.Date.of.IP.Administration` <- 
  convert_excel_date(br_df$`Start.Date.of.IP.Administration`)
br_df$`date.of.progression` <- 
  convert_excel_date(br_df$`date.of.progression`)
br_df$`date.of.death.or.last.FU` <- 
  convert_excel_date(br_df$`date.of.death.or.last.FU`)
# br_df$`Basline.Biopsy` <- convert_excel_date(br_df$`Basline.Biopsy`)
```

``` r
br_df$`Start.Date.of.IP.Administration`[1:5]
```

    ## [1] "2024-02-28" "2024-05-07" "2024-05-22" "2024-05-29" "2024-06-04"

``` r
br_df[c("Subject.Enrollment.Number", "Start.Date.of.IP.Administration", "date.of.progression", "date.of.death.or.last.FU")]
```

    ##    Subject.Enrollment.Number Start.Date.of.IP.Administration
    ## 1                 B11-A-1-01                      2024-02-28
    ## 4                 B11-A-1-04                      2024-05-07
    ## 5                 B11-A-1-05                      2024-05-22
    ## 6                 B11-A-1-06                      2024-05-29
    ## 7                 B11-A-1-07                      2024-06-04
    ## 9                 B11-A-1-09                      2024-06-12
    ## 10                B11-A-1-10                      2024-06-12
    ## 11                B11-A-1-11                      2024-06-26
    ## 12                B11-A-1-12                      2024-07-02
    ## 13                B11-A-1-13                      2024-06-25
    ## 14                B11-A-1-14                      2024-07-17
    ## 17                B11-A-1-17                      2024-07-16
    ## 18                B11-A-1-18                      2024-07-24
    ## 19                B11-A-1-19                      2024-09-10
    ## 20                B11-A-1-20                      2024-11-13
    ## 21                B11-A-1-21                      2024-11-27
    ## 22                B11-A-1-22                      2024-11-28
    ## 23                B11-A-1-23                      2024-11-27
    ## 24                B11-A-1-24                      2024-12-11
    ## 26                B11-A-1-26                      2024-12-19
    ## 27                B11-A-1-27                      2024-12-24
    ## 28                B11-A-1-28                      2025-01-15
    ## 29                B11-A-1-29                      2025-02-04
    ## 30                B11-A-1-30                      2025-02-11
    ## 31                B11-A-1-31                      2025-02-12
    ## 32                B11-A-1-32                      2025-02-03
    ## 33                B11-A-1-33                      2025-02-12
    ## 34                B11-A-1-34                      2025-03-25
    ##    date.of.progression date.of.death.or.last.FU
    ## 1           2024-04-25               2024-08-17
    ## 4           2024-06-12               2024-07-05
    ## 5           2024-07-17               2024-08-15
    ## 6           2024-07-03               2024-07-28
    ## 7           2025-07-23               2025-07-23
    ## 9           2025-07-23               2025-07-23
    ## 10          2024-10-01               2025-01-06
    ## 11          2025-01-27               2025-02-17
    ## 12          2025-07-23               2025-07-23
    ## 13          2024-08-06               2024-09-08
    ## 14          2024-09-25               2025-01-09
    ## 17          2024-09-04               2025-01-08
    ## 18          2025-07-23               2025-07-23
    ## 19          2025-07-23               2025-07-23
    ## 20          2025-02-25               2025-07-23
    ## 21          2025-07-23               2025-07-23
    ## 22          2025-01-23               2025-07-23
    ## 23          2025-05-12               2025-07-23
    ## 24          2025-03-12               2025-04-12
    ## 26          2025-05-08               2025-07-23
    ## 27          2025-02-19               2025-05-19
    ## 28          2025-03-17               2025-07-23
    ## 29          2025-07-23               2025-07-23
    ## 30          2025-07-23               2025-07-23
    ## 31          2025-03-27               2025-04-02
    ## 32          2025-07-23               2025-07-23
    ## 33          2025-07-23               2025-07-23
    ## 34          2025-07-23               2025-07-23

If the `date of progression` is not available, we will use the
`date of death or last follow-up` to calculate the PFS.

``` r
# Create the PFS and OS columns by subtracting the Baseline Biopsy date from the `date.of.progression` and `date.of.death.or.last.FU`
br_df$OS <- as.numeric(difftime(br_df$`date.of.death.or.last.FU`, br_df$`Start.Date.of.IP.Administration`, units = "days"))
br_df$OS
```

    ##  [1] 171  59  85  60 414 406 208 236 386  75 176 176 364 316 252 238 237 238 122
    ## [20] 216 146 189 169 162  49 170 161 120

``` r
# Count the number of NA values in the OS column
sum(is.na(br_df$OS))  # Count the number of NA values in the OS
```

    ## [1] 0

``` r
# Calculate PFS values
for (i in seq_len(nrow(br_df))) {
  prog_date <- br_df$`date.of.progression`[i]
  start_date <- br_df$`Start.Date.of.IP.Administration`[i]
  fu_date <- br_df$`date.of.death.or.last.FU`[i]
  
  if (is.na(prog_date)) {
    # Use follow-up date if progression date is missing
    br_df$PFS[i] <- as.numeric(difftime(fu_date, start_date, 
                                         units = "days"))
  } else {
    # Use progression date if available
    br_df$PFS[i] <- as.numeric(difftime(prog_date, start_date, 
                                         units = "days"))
  }
}

br_df$PFS
```

    ##  [1]  57  36  56  35 414 406 111 215 386  42  70  50 364 316 104 238  56 166  91
    ## [20] 140  57  61 169 162  43 170 161 120

``` r
nrow(br_df)
```

    ## [1] 28

``` r
# Count the number of NA values in the OS or PFS columns
sum(is.na(br_df$OS) | is.na(br_df$PFS))  # Count the number of NA values in the OS or PFS columns
```

    ## [1] 0

``` r
# Print "Subject.Enrollment.Number", "PFS", and "OS" columns
br_df[c("Subject.Enrollment.Number", "Start.Date.of.IP.Administration", "date.of.progression", "PFS", "date.of.death.or.last.FU", "OS")]
```

    ##    Subject.Enrollment.Number Start.Date.of.IP.Administration
    ## 1                 B11-A-1-01                      2024-02-28
    ## 4                 B11-A-1-04                      2024-05-07
    ## 5                 B11-A-1-05                      2024-05-22
    ## 6                 B11-A-1-06                      2024-05-29
    ## 7                 B11-A-1-07                      2024-06-04
    ## 9                 B11-A-1-09                      2024-06-12
    ## 10                B11-A-1-10                      2024-06-12
    ## 11                B11-A-1-11                      2024-06-26
    ## 12                B11-A-1-12                      2024-07-02
    ## 13                B11-A-1-13                      2024-06-25
    ## 14                B11-A-1-14                      2024-07-17
    ## 17                B11-A-1-17                      2024-07-16
    ## 18                B11-A-1-18                      2024-07-24
    ## 19                B11-A-1-19                      2024-09-10
    ## 20                B11-A-1-20                      2024-11-13
    ## 21                B11-A-1-21                      2024-11-27
    ## 22                B11-A-1-22                      2024-11-28
    ## 23                B11-A-1-23                      2024-11-27
    ## 24                B11-A-1-24                      2024-12-11
    ## 26                B11-A-1-26                      2024-12-19
    ## 27                B11-A-1-27                      2024-12-24
    ## 28                B11-A-1-28                      2025-01-15
    ## 29                B11-A-1-29                      2025-02-04
    ## 30                B11-A-1-30                      2025-02-11
    ## 31                B11-A-1-31                      2025-02-12
    ## 32                B11-A-1-32                      2025-02-03
    ## 33                B11-A-1-33                      2025-02-12
    ## 34                B11-A-1-34                      2025-03-25
    ##    date.of.progression PFS date.of.death.or.last.FU  OS
    ## 1           2024-04-25  57               2024-08-17 171
    ## 4           2024-06-12  36               2024-07-05  59
    ## 5           2024-07-17  56               2024-08-15  85
    ## 6           2024-07-03  35               2024-07-28  60
    ## 7           2025-07-23 414               2025-07-23 414
    ## 9           2025-07-23 406               2025-07-23 406
    ## 10          2024-10-01 111               2025-01-06 208
    ## 11          2025-01-27 215               2025-02-17 236
    ## 12          2025-07-23 386               2025-07-23 386
    ## 13          2024-08-06  42               2024-09-08  75
    ## 14          2024-09-25  70               2025-01-09 176
    ## 17          2024-09-04  50               2025-01-08 176
    ## 18          2025-07-23 364               2025-07-23 364
    ## 19          2025-07-23 316               2025-07-23 316
    ## 20          2025-02-25 104               2025-07-23 252
    ## 21          2025-07-23 238               2025-07-23 238
    ## 22          2025-01-23  56               2025-07-23 237
    ## 23          2025-05-12 166               2025-07-23 238
    ## 24          2025-03-12  91               2025-04-12 122
    ## 26          2025-05-08 140               2025-07-23 216
    ## 27          2025-02-19  57               2025-05-19 146
    ## 28          2025-03-17  61               2025-07-23 189
    ## 29          2025-07-23 169               2025-07-23 169
    ## 30          2025-07-23 162               2025-07-23 162
    ## 31          2025-03-27  43               2025-04-02  49
    ## 32          2025-07-23 170               2025-07-23 170
    ## 33          2025-07-23 161               2025-07-23 161
    ## 34          2025-07-23 120               2025-07-23 120

``` r
# Drop rows with NA in PFS or OS
br_df <- br_df[!is.na(br_df$PFS) & !is.na(br_df$OS), ]
nrow(br_df)  # Check the number of rows after dropping NA values
```

    ## [1] 28

``` r
# Drop rows with negative PFS or OS values
# br_df <- br_df[br_df$PFS >= 0 & br_df$OS >= 0, ]
# nrow(br_df)  # Check the number of rows after dropping negative PFS or OS values
```

``` r
# Convert the PFS and OS columns to months
br_df$PFS <- br_df$PFS / 30.44  # Convert PFS from days to months
br_df$OS <- br_df$OS / 30.44  # Convert OS from days to months
```

``` r
br_df$`SNV/Indel`
```

    ##  [1] "PIK3CA"                 "KRAS G12D"              NA                      
    ##  [4] "KRAS G12D\nPIK3CA "     "ATM "                   "PIK3CA\nARID1A    "    
    ##  [7] NA                       NA                       "PIK3CA\nARID1A\nTP53 " 
    ## [10] "TSC2 "                  "TSC2 "                  NA                      
    ## [13] "TP53"                   NA                       NA                      
    ## [16] "PIK3CA  "               NA                       "TP53\nARID1A"          
    ## [19] NA                       "ATM\nARID1A"            NA                      
    ## [22] "ATM del"                "KRAS\nPIK3CA\nPTEN del" "PIK3CA"                
    ## [25] NA                       "MAP2K1\nTSC1 del"       NA                      
    ## [28] "ATM del"

A function to refine the `SNV/Indel` column by grouping values by user
input.

``` r
library(stringr)
refine_snv_indel <- function(br_df) {
  # Load required package
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required. Please install it.")
  }

  # Copy the original data
  df <- br_df
  # Backup the original column for comparison
  original_values <- df$`SNV/Indel`
  
  # Create a working column to modify
  df$`SNV/Indel` <- as.character(df$`SNV/Indel`)  # Ensure character type
  converted_idx <- integer(0)  # Store indices of converted rows
  
  repeat {
    # Get unconverted (and non-NA) values
    remaining_idx <- setdiff(which(!is.na(df$`SNV/Indel`)), converted_idx)
    
    if (length(remaining_idx) == 0) {
      cat("All values have been converted.\n")
      break
    }
    
    # Show current remaining values to the user
    cat("\nRemaining `SNV/Indel` values:\n")
    print(unique(df$`SNV/Indel`[remaining_idx]))
    
    # Prompt user input
    user_input <- readline(prompt = "\nEnter string to standardize (Enter to finish): ")
    
    # Exit if user just presses enter
    if (user_input == "") {
      cat("No input detected. Exiting refinement.\n")
      break
    }
    
    # Find matching indices
    matched_idx <- remaining_idx[
      stringr::str_detect(df$`SNV/Indel`[remaining_idx], fixed(user_input))
    ]
    
    # If nothing matched, inform user and continue
    if (length(matched_idx) == 0) {
      cat("No values matched your input.\n")
      next
    }
    
    # Replace matched values with the user input
    df$`SNV/Indel`[matched_idx] <- user_input
    
    # Update converted indices
    converted_idx <- union(converted_idx, matched_idx)
    
    # Summary of replacements
    cat(sprintf("Replaced %d entries with \"%s\"\n", length(matched_idx), user_input))
    
    # Optional early exit if all rows are converted
    if (setequal(converted_idx, which(!is.na(df$`SNV/Indel`)))) {
      cat("All non-NA entries have been grouped by user input.\n")
      break
    }
  }

  # Return a summary table with before and after values
  result_df <- br_df
  result_df$SNV_Indel_before <- original_values
  result_df$SNV_Indel_after <- df$`SNV/Indel`
  
  return(result_df)
}
```

``` r
# Apply the refine_snv_indel function to br_df
df_snv_grouped <- refine_snv_indel(br_df)
```

    ## 
    ## Remaining `SNV/Indel` values:
    ##  [1] "PIK3CA"                 "KRAS G12D"              "KRAS G12D\nPIK3CA "    
    ##  [4] "ATM "                   "PIK3CA\nARID1A    "     "PIK3CA\nARID1A\nTP53 " 
    ##  [7] "TSC2 "                  "TP53"                   "PIK3CA  "              
    ## [10] "TP53\nARID1A"           "ATM\nARID1A"            "ATM del"               
    ## [13] "KRAS\nPIK3CA\nPTEN del" "MAP2K1\nTSC1 del"      
    ## 
    ## Enter string to standardize (Enter to finish): 
    ## No input detected. Exiting refinement.

``` r
df_snv_grouped[c("SNV_Indel_before", "SNV_Indel_after")]
```

    ##          SNV_Indel_before        SNV_Indel_after
    ## 1                  PIK3CA                 PIK3CA
    ## 4               KRAS G12D              KRAS G12D
    ## 5                    <NA>                   <NA>
    ## 6      KRAS G12D\nPIK3CA      KRAS G12D\nPIK3CA 
    ## 7                    ATM                    ATM 
    ## 9      PIK3CA\nARID1A         PIK3CA\nARID1A    
    ## 10                   <NA>                   <NA>
    ## 11                   <NA>                   <NA>
    ## 12  PIK3CA\nARID1A\nTP53   PIK3CA\nARID1A\nTP53 
    ## 13                  TSC2                   TSC2 
    ## 14                  TSC2                   TSC2 
    ## 17                   <NA>                   <NA>
    ## 18                   TP53                   TP53
    ## 19                   <NA>                   <NA>
    ## 20                   <NA>                   <NA>
    ## 21               PIK3CA                 PIK3CA  
    ## 22                   <NA>                   <NA>
    ## 23           TP53\nARID1A           TP53\nARID1A
    ## 24                   <NA>                   <NA>
    ## 26            ATM\nARID1A            ATM\nARID1A
    ## 27                   <NA>                   <NA>
    ## 28                ATM del                ATM del
    ## 29 KRAS\nPIK3CA\nPTEN del KRAS\nPIK3CA\nPTEN del
    ## 30                 PIK3CA                 PIK3CA
    ## 31                   <NA>                   <NA>
    ## 32       MAP2K1\nTSC1 del       MAP2K1\nTSC1 del
    ## 33                   <NA>                   <NA>
    ## 34                ATM del                ATM del

``` r
# Remove NA values in the SNV/Indel column
# df_snv_grouped <- df_snv_grouped[!is.na(df_snv_grouped$SNV_Indel_after), ]
# df_snv_grouped[c("SNV_Indel_before", "SNV_Indel_after")]
```

# Copy the modified df

``` r
# Copy the modified data frame to df
df <- df_snv_grouped
# Prepare an list of columns to keep for survival analysis
df_surv_cols <- c("PFS", "PD", "OS", "Death")
```

``` r
# Drop the original SNV/Indel column
df$`SNV/Indel` <- NULL
df$`SNV_Indel_before` <- NULL
# Append "SNV_Indel_after" to df_surv_cols
df_surv_cols <- c(df_surv_cols, "SNV_Indel_after")
# Show columns in df not in df_surv_cols
setdiff(colnames(df), df_surv_cols)
```

    ##  [1] "No."                                  
    ##  [2] "Subject.Enrollment.Number"            
    ##  [3] "Best.Response.(Cycle)"                
    ##  [4] "HER2"                                 
    ##  [5] "EBV"                                  
    ##  [6] "MLH"                                  
    ##  [7] "PD-L1(28-8)"                          
    ##  [8] "PD-L1(22C3)"                          
    ##  [9] "CLDN18"                               
    ## [10] "TMB"                                  
    ## [11] "MSI"                                  
    ## [12] "MMRD"                                 
    ## [13] "CNV"                                  
    ## [14] "Fusion"                               
    ## [15] "age"                                  
    ## [16] "sex"                                  
    ## [17] "status"                               
    ## [18] "reason.for.end.of.administration"     
    ## [19] "Subject.Enrollment.Number.1"          
    ## [20] "Subject.Initial"                      
    ## [21] "Base.line.TL.Sum.(cm)"                
    ## [22] "Max.change.from.baseline.(Percentage)"
    ## [23] "Last.Date.of.IP.Administration"       
    ## [24] "Durvalumab.(fix.dose)"                
    ## [25] "First.dose.(Ceralasertib)"            
    ## [26] "Current.Dose.(Ceralasertib)"          
    ## [27] "DR"                                   
    ## [28] "DR.Cycle"                             
    ## [29] "DR.date"                              
    ## [30] "Reason.of.Dose.reduction"             
    ## [31] "Start.Date.of.IP.Administration"      
    ## [32] "date.of.progression"                  
    ## [33] "date.of.death.or.last.FU"             
    ## [34] "Basline.Biopsy"                       
    ## [35] "Post-Biopsy.(C1D8.pre)"               
    ## [36] "Post-Biopsy.(C2D1.pre)"               
    ## [37] "Basline"                              
    ## [38] "Post-Biopsy.(C1D8.pre).1"             
    ## [39] "Post-Biopsy.(C2D1.pre).1"             
    ## [40] "BL"                                   
    ## [41] "FU1"                                  
    ## [42] "FU2"                                  
    ## [43] "BL.1"                                 
    ## [44] "FU1.1"                                
    ## [45] "FU2.1"                                
    ## [46] "FU3"                                  
    ## [47] "FU4"                                  
    ## [48] "FU5"                                  
    ## [49] "FU6"                                  
    ## [50] "FU7"                                  
    ## [51] "FU8"                                  
    ## [52] "FU9"                                  
    ## [53] "FU10"                                 
    ## [54] "FU11"                                 
    ## [55] "FU12"                                 
    ## [56] "FU13"                                 
    ## [57] "FU14"

## Preprocess other columns.

``` r
table(df$HER2)
```

    ## 
    ##        0       1+ negative 
    ##       15       12        1

``` r
df$HER2 <- ifelse(df$HER2 == "negative", "negative", "neutral")
df$HER2 <- factor(df$HER2, unique(df$HER2[!is.na(df$HER2)]))
table(df$HER2)  # Check the distribution of HER2 after preprocessing
```

    ## 
    ##  neutral negative 
    ##       27        1

``` r
df_surv_cols <- c(df_surv_cols, "HER2")  # Append HER2 to df_surv_cols
```

``` r
table(df$EBV)
```

    ## 
    ##            neg Not applicable            pos 
    ##             12              1              8

``` r
df$EBV <- factor(df$EBV, unique(df$EBV[!is.na(df$EBV)]))
table(df$EBV)  # Check the distribution of EBV after preprocessing
```

    ## 
    ##            pos            neg Not applicable 
    ##              8             12              1

``` r
df_surv_cols <- c(df_surv_cols, "EBV")  # Append EBV to df_surv_cols
```

``` r
table(df$MLH)
```

    ## 
    ##                                     Intact 
    ##                                         22 
    ## Loss of expression in 100% of tumor cells  
    ##                                          4 
    ##  Loss of expression in 95% of tumor cells  
    ##                                          1

``` r
# Concatenate 'Loss of expression in 95% of tumor cells' and 'Loss of expression in 100% of tumor cells' into 'Loss of expression in over 95% of tumor cells'
df$MLH <- ifelse(df$MLH != "Intact", "Loss of expression in over 95% of tumor cells", df$MLH)
table(df$MLH)  # Check the distribution of MLH after preprocessing
```

    ## 
    ##                                        Intact 
    ##                                            22 
    ## Loss of expression in over 95% of tumor cells 
    ##                                             5

``` r
df_surv_cols <- c(df_surv_cols, "MLH")  # Append MLH to df_surv_cols
```

``` r
table(df$`PD-L1(28-8)`)
```

    ## 
    ##  1  2  3  5 10 35 60 
    ##  2  1  6  7  3  1  1

``` r
typeof(df$`PD-L1(28-8)`)  # Check the type of the PD-L1(28-8) column
```

    ## [1] "double"

``` r
df_surv_cols <- c(df_surv_cols, "PD-L1(28-8)")  # Append PD-L1(28-8) to df_surv_cols
df_surv_cols
```

    ## [1] "PFS"             "PD"              "OS"              "Death"          
    ## [5] "SNV_Indel_after" "HER2"            "EBV"             "MLH"            
    ## [9] "PD-L1(28-8)"

``` r
table(df$`PD-L1(22C3)`)
```

    ## 
    ##  2  3  5 10 20 25 90 
    ##  1  3  8  3  2  1  1

``` r
typeof(df$`PD-L1(22C3)`)  # Check the type of the PD-L1(22C3) column
```

    ## [1] "double"

``` r
df_surv_cols <- c(df_surv_cols, "PD-L1(22C3)")  # Append PD-L1(22C3) to df_surv_cols
df_surv_cols
```

    ##  [1] "PFS"             "PD"              "OS"              "Death"          
    ##  [5] "SNV_Indel_after" "HER2"            "EBV"             "MLH"            
    ##  [9] "PD-L1(28-8)"     "PD-L1(22C3)"

``` r
table(df$TMB)
```

    ## 
    ## High  Low 
    ##    8   13

``` r
df$TMB <- factor(df$TMB, levels = c("High", "Low"))
table(df$TMB)  # Check the distribution of TMB after preprocessing
```

    ## 
    ## High  Low 
    ##    8   13

``` r
df_surv_cols <- c(df_surv_cols, "TMB")  # Append TMB to df_surv_cols
df_surv_cols
```

    ##  [1] "PFS"             "PD"              "OS"              "Death"          
    ##  [5] "SNV_Indel_after" "HER2"            "EBV"             "MLH"            
    ##  [9] "PD-L1(28-8)"     "PD-L1(22C3)"     "TMB"

``` r
table(df$MSI)
```

    ## 
    ##          High Indeterminate           MSS 
    ##             4             1            18

``` r
df$MSI <- factor(df$MSI, levels = unique(df$MSI[!is.na(df$MSI)]))
table(df$MSI)  # Check the distribution of MSI after preprocessing
```

    ## 
    ##           MSS Indeterminate          High 
    ##            18             1             4

``` r
df_surv_cols <- c(df_surv_cols, "MSI")  # Append MSI to df_surv_cols
df_surv_cols
```

    ##  [1] "PFS"             "PD"              "OS"              "Death"          
    ##  [5] "SNV_Indel_after" "HER2"            "EBV"             "MLH"            
    ##  [9] "PD-L1(28-8)"     "PD-L1(22C3)"     "TMB"             "MSI"

``` r
table(df$CNV)
```

    ## 
    ## CCND3 amp\nVEGFA amp            CCNE1 amp            CDK6  amp 
    ##                    1                    1                    1 
    ##            FGFR2 amp             MDM2 amp              MET amp 
    ##                    1                    2                    1 
    ##           PIK3CA amp           RICTOR amp            SMAD3 Del 
    ##                    1                    2                    1

``` r
df$CNV <- factor(df$CNV, levels = unique(df$CNV[!is.na(df$CNV)]))
table(df$CNV)  # Check the distribution of CNV after preprocessing
```

    ## 
    ##           RICTOR amp           PIK3CA amp            CDK6  amp 
    ##                    2                    1                    1 
    ## CCND3 amp\nVEGFA amp            CCNE1 amp             MDM2 amp 
    ##                    1                    1                    2 
    ##            FGFR2 amp            SMAD3 Del              MET amp 
    ##                    1                    1                    1

``` r
# Drop the df$CNV column if it is not needed
df$CNV <- NULL
```

``` r
table(df$sex)
```

    ## 
    ##  F  M 
    ##  7 21

``` r
df$sex <- factor(df$sex, levels = unique(df$sex[!is.na(df$sex)]))
table(df$sex)  # Check the distribution of sex after preprocessing
```

    ## 
    ##  M  F 
    ## 21  7

``` r
df$MMRD <- factor(df$MMRD, levels = unique(df$MMRD[!is.na(df$MMRD)]))
table(df$MMRD)
```

    ## 
    ## MMRp MMRd 
    ##   22    6

``` r
df_surv_cols <- c(df_surv_cols, "MMRD")  # Append MSI to df_surv_cols
df_surv_cols
```

    ##  [1] "PFS"             "PD"              "OS"              "Death"          
    ##  [5] "SNV_Indel_after" "HER2"            "EBV"             "MLH"            
    ##  [9] "PD-L1(28-8)"     "PD-L1(22C3)"     "TMB"             "MSI"            
    ## [13] "MMRD"

``` r
table(df$`Best.Response.(Cycle)`)
```

    ## 
    ##      PD PR (C3) PR (C5) PR (C6)  PR(C3)      SD 
    ##      12       5       2       1       1       7

``` r
df_surv_cols <- c(df_surv_cols, 'Best.Response.(Cycle)')
df_surv_cols
```

    ##  [1] "PFS"                   "PD"                    "OS"                   
    ##  [4] "Death"                 "SNV_Indel_after"       "HER2"                 
    ##  [7] "EBV"                   "MLH"                   "PD-L1(28-8)"          
    ## [10] "PD-L1(22C3)"           "TMB"                   "MSI"                  
    ## [13] "MMRD"                  "Best.Response.(Cycle)"

``` r
df_surv_cols <- c(df_surv_cols, "sex", "age")
df_surv_cols
```

    ##  [1] "PFS"                   "PD"                    "OS"                   
    ##  [4] "Death"                 "SNV_Indel_after"       "HER2"                 
    ##  [7] "EBV"                   "MLH"                   "PD-L1(28-8)"          
    ## [10] "PD-L1(22C3)"           "TMB"                   "MSI"                  
    ## [13] "MMRD"                  "Best.Response.(Cycle)" "sex"                  
    ## [16] "age"

# Factorize df into two levels: ATM and non-ATM.

## Show the SNV_Indel-groups

``` r
# Check the dimensions of the data frame after preprocessing
# library(skimr)
# skim(df_surv)  # Get a summary of the data frame

table(df$SNV_Indel_after, useNA = 'ifany')
```

    ## 
    ##                   ATM                 ATM del            ATM\nARID1A 
    ##                      1                      2                      1 
    ##              KRAS G12D     KRAS G12D\nPIK3CA  KRAS\nPIK3CA\nPTEN del 
    ##                      1                      1                      1 
    ##       MAP2K1\nTSC1 del                 PIK3CA               PIK3CA   
    ##                      1                      2                      1 
    ##     PIK3CA\nARID1A      PIK3CA\nARID1A\nTP53                    TP53 
    ##                      1                      1                      1 
    ##           TP53\nARID1A                  TSC2                    <NA> 
    ##                      1                      2                     11

``` r
nrow(df)
```

    ## [1] 28

## Preprocess the SNV_Indel_after

``` r
# Keep only the relevant columns for survival analysis
df_surv <- df[, df_surv_cols]

# Exclude certain columns but preserve NA
MAKE_NA <- c()

df_surv <- df_surv %>%
  filter(!SNV_Indel_after %in% MAKE_NA)
table(df_surv$SNV_Indel_after, useNA = 'ifany')  # Check the distribution of SNV_Indel_after)
```

    ## 
    ##                   ATM                 ATM del            ATM\nARID1A 
    ##                      1                      2                      1 
    ##              KRAS G12D     KRAS G12D\nPIK3CA  KRAS\nPIK3CA\nPTEN del 
    ##                      1                      1                      1 
    ##       MAP2K1\nTSC1 del                 PIK3CA               PIK3CA   
    ##                      1                      2                      1 
    ##     PIK3CA\nARID1A      PIK3CA\nARID1A\nTP53                    TP53 
    ##                      1                      1                      1 
    ##           TP53\nARID1A                  TSC2                    <NA> 
    ##                      1                      2                     11

``` r
nrow(df_surv)  # Check the number of rows after filtering
```

    ## [1] 28

``` r
# Divide values into two groups - ATM_mut+KRAS and otherwise 
# in a new column 'ATM_mut' based on the 'SNV_Indel_after' column
df_surv <- df_surv %>%
  mutate(ATM_mut = ifelse(is.na(SNV_Indel_after), "Otherwise",
                                    ifelse(SNV_Indel_after == "ATM" |
                                             SNV_Indel_after == "ARID1A" |
                                             SNV_Indel_after == "TP53",
                                     "ATM_mut_included", "Otherwise"))) %>%
  group_by(ATM_mut)


# If only want to use ATM only
# df_surv <- df_surv %>%
#   mutate(ATM_mut = ifelse(is.na(SNV_Indel_after), "Otherwise", 
#                                     ifelse(SNV_Indel_after == "ATM",
#                                      "ATM_mut_included", "Otherwise"))) %>%
#   group_by(ATM_mut)



df_surv$SNV_Indel_after <- NULL  # Remove the SNV_Indel_after column as it is not needed for survival analysis

table(df_surv$ATM_mut)  # Check the distribution of ATM_mut
```

    ## 
    ## ATM_mut_included        Otherwise 
    ##                1               27

## PFS

``` r
km_fit_group <- survfit(Surv(PFS, PD) ~ ATM_mut, data = df_surv)

# Solve the difference between the mean time and number of samples
# Calculate the mean time for each group
mean_time <- mean(df_surv$PFS[df_surv$ATM_mut == "ATM_mut_included"], na.rm = TRUE)
mean_time2 <- mean(df_surv$PFS[df_surv$ATM_mut == "Otherwise"], na.rm = TRUE)
# Calculate the number of samples for each group
n_group1 <- nrow(df_surv[df_surv$ATM_mut == "ATM_mut_included", ])
n_group2 <- nrow(df_surv[df_surv$ATM_mut == "Otherwise", ])
n_group_text <- paste0("ATM_mut_included: ", n_group1, "\nOtherwise: ", n_group2)

mean_time <- round(mean_time, 2)
mean_time2 <- round(mean_time2, 2)
diff_time <- mean_time - mean_time2
# Create a text label for the mean time
mean_time_text <- paste0("mPFS of ATMmut: ", mean_time, " months\n",
                          "mPFS of non-ATMmut: ", mean_time2, " months")
# Create a text label for the difference in mean time
diff_time_text <- paste0("Difference in mean time: ", diff_time, " months")
# Create a text label for the number of samples

n_group_text <- paste0("ATMmut: ", n_group1, "\nnon-ATMmut: ", n_group2)
medians <- summary(km_fit_group)$table[, "median"]


# Calculate the log-rank test and p-value
logrank_test <- survdiff(Surv(PFS, PD) ~ ATM_mut, data = df_surv)
pval <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
pval_text <- paste0("p = ", format.pval(pval, digits = 3, eps = .001))

# Total number of samples
total_n <- nrow(df_surv)
n_text <- paste0("n = ", total_n)


plot_obj <- ggsurvplot(
  km_fit_group,
  data = df_surv,
  xlab = "Time (months)",
  ylab = "Survival Probability",
  title = "PFS of ATMmut vs. non-ATMmut",
  palette = c("blue", "red"),
  conf.int = FALSE,
  pval = pval_text,
  pval.coord = c(10, 0.05),
  legend.title = "Group",
  legend.labs = c("ATMmut", "nonATMmut"),
  risk.table = TRUE,
  risk.table.title = NULL,
  risk.table.y.text = FALSE,
  break.time.by = 10,
  xlim = c(0, 40),
  surv.scale = "percent",
  mean_time_text = mean_time_text,
  diff_time_text = diff_time_text,
  ggtheme = theme_minimal()
)

# Add the total number of samples to the plot
plot_obj$plot <- plot_obj$plot +
  annotate("text", x = 35, y = 1.0, label = n_text, hjust = 1, size = 4) +
  scale_x_continuous(breaks = seq(0, 40, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = scales::percent_format(accuracy = 1)) +
  annotate("text", x = 30, y = 0.8, label = mean_time_text, hjust = 1, size = 4) +
  annotate("text", x = 30, y = 0.55, label = diff_time_text, hjust = 1, size = 4) +
  # Dotted Line for 50% survival
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "#306c81", linewidth = 0.5) +
  geom_vline(xintercept = medians[1], linetype = "dashed", color = "blue", size = 0.5) +
  geom_vline(xintercept = medians[2], linetype = "dashed", color = "red", size = 0.5) +
  annotate("text", x = medians[1], y = 0.1, label = paste0("Med.: ", round(medians[1], 2)), color = "blue", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = medians[2], y = 0.1, label = paste0("Med.: ", round(medians[2], 2)), color = "red", angle = 90, vjust = -0.5, size = 3)
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

### Plot

``` r
# Print the plot
print(plot_obj)
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_vline()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_vline()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text()`).

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-64-1.png)<!-- -->

``` r
# Save the plot into 1024*768 png file to OUTPUT_DIR
# png(file.path(OUTPUT_DIR, "[VT]PFS_surv_grouped_by_ATMmut.png"), width = 2000, height = 1500, res = 300)
# 
# plot_obj
# 
# dev.off()
```

## OS

``` r
km_fit_group <- survfit(Surv(OS, Death) ~ ATM_mut, data = df_surv)

# Solve the difference between the mean time and number of samples
# Calculate the mean time for each group
mean_time <- mean(df_surv$OS[df_surv$ATM_mut == "ATM_mut_included"], na.rm = TRUE)
mean_time2 <- mean(df_surv$OS[df_surv$ATM_mut == "Otherwise"], na.rm = TRUE)
# Calculate the number of samples for each group
n_group1 <- nrow(df_surv[df_surv$ATM_mut == "ATM_mut_included", ])
n_group2 <- nrow(df_surv[df_surv$ATM_mut == "Otherwise", ])
n_group_text <- paste0("ATM_mut_included: ", n_group1, "\nOtherwise: ", n_group2)

mean_time <- round(mean_time, 2)
mean_time2 <- round(mean_time2, 2)
diff_time <- mean_time - mean_time2
# Create a text label for the mean time
mean_time_text <- paste0("mOS of ATMmut: ", mean_time, " months\n",
                          "mOS of non-ATMmut: ", mean_time2, " months")
# Create a text label for the difference in mean time
diff_time_text <- paste0("Difference in mean time: ", diff_time, " months")
# Create a text label for the number of samples

n_group_text <- paste0("ATMmut: ", n_group1, "\nnon-ATMmut: ", n_group2)
medians <- summary(km_fit_group)$table[, "median"]


# Calculate the log-rank test and p-value
logrank_test <- survdiff(Surv(OS, Death) ~ ATM_mut, data = df_surv)
pval <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
pval_text <- paste0("p = ", format.pval(pval, digits = 3, eps = .001))

# Total number of samples
total_n <- nrow(df_surv)
n_text <- paste0("n = ", total_n)


plot_obj <- ggsurvplot(
  km_fit_group,
  data = df_surv,
  xlab = "Time (months)",
  ylab = "Survival Probability",
  title = "OS of ATMmut vs. non-ATMmut",
  palette = c("blue", "red"),
  conf.int = FALSE,
  pval = pval_text,
  pval.coord = c(10, 0.05),
  legend.title = "Group",
  legend.labs = c("ATMmut", "non-ATMmut"),
  risk.table = TRUE,
  risk.table.title = NULL,
  risk.table.y.text = FALSE,
  break.time.by = 10,
  xlim = c(0, 40),
  surv.scale = "percent",
  mean_time_text = mean_time_text,
  diff_time_text = diff_time_text,
  ggtheme = theme_minimal()
)

# Add the total number of samples to the plot
plot_obj$plot <- plot_obj$plot +
  annotate("text", x = 35, y = 1.0, label = n_text, hjust = 1, size = 4) +
  scale_x_continuous(breaks = seq(0, 40, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = scales::percent_format(accuracy = 1)) +
  annotate("text", x = 30, y = 0.8, label = mean_time_text, hjust = 1, size = 4) +
  annotate("text", x = 30, y = 0.55, label = diff_time_text, hjust = 1, size = 4) +
  # Dotted Line for 50% survival
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "#306c81", linewidth = 0.5) +
  geom_vline(xintercept = medians[1], linetype = "dashed", color = "blue", size = 0.5) +
  geom_vline(xintercept = medians[2], linetype = "dashed", color = "red", size = 0.5) +
  annotate("text", x = medians[1], y = 0.1, label = paste0("Med.: ", round(medians[1], 2)), color = "blue", angle = 90, vjust = -0.5, size = 3) +
  annotate("text", x = medians[2], y = 0.1, label = paste0("Med.: ", round(medians[2], 2)), color = "red", angle = 90, vjust = -0.5, size = 3)
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

### Plot

``` r
# Print the plot
print(plot_obj)
```

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_vline()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_vline()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_text()`).

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

``` r
# Save the plot into 1024*768 png file to OUTPUT_DIR
# png(file.path(OUTPUT_DIR, "[VT]OS_surv_grouped_by_ATMmut.png"), width = 2000, height = 1500, res = 300)
# 
# plot_obj
# 
# dev.off()
```

## Cox analysis

``` r
head(df_surv)
```

    ## # A tibble: 6 × 16
    ## # Groups:   ATM_mut [1]
    ##     PFS    PD    OS Death HER2     EBV   MLH   `PD-L1(28-8)` `PD-L1(22C3)` TMB  
    ##   <dbl> <dbl> <dbl> <dbl> <fct>    <fct> <chr>         <dbl>         <dbl> <fct>
    ## 1  1.87     1  5.62     1 neutral  pos   Inta…             5            20 Low  
    ## 2  1.18     1  1.94     1 negative <NA>  Inta…             1            NA Low  
    ## 3  1.84     1  2.79     1 neutral  neg   Inta…             3             5 Low  
    ## 4  1.15     1  1.97     1 neutral  neg   Inta…            35             5 High 
    ## 5 13.6      0 13.6      0 neutral  neg   Inta…             1            NA Low  
    ## 6 13.3      0 13.3      0 neutral  pos   Inta…            10            20 Low  
    ## # ℹ 6 more variables: MSI <fct>, MMRD <fct>, `Best.Response.(Cycle)` <chr>,
    ## #   sex <fct>, age <dbl>, ATM_mut <chr>

Create new columns: Responder(SD, PD are Non-responder), DCB(Durable
Clinical Benefit)

``` r
df_surv <- df_surv |> mutate(Responder = ifelse(`Best.Response.(Cycle)` == "SD" |
                                           `Best.Response.(Cycle)` == "PD", 0, 1),
                             DCB = ifelse(grepl("PR", `Best.Response.(Cycle)`) |
                                            (`Best.Response.(Cycle)` == "SD" & PFS >= 6),
                                          1, 0))
table(df_surv$Responder)
```

    ## 
    ##  0  1 
    ## 19  9

``` r
table(df_surv$DCB)
```

    ## 
    ##  0  1 
    ## 17 11

``` r
colnames(df_surv)
```

    ##  [1] "PFS"                   "PD"                    "OS"                   
    ##  [4] "Death"                 "HER2"                  "EBV"                  
    ##  [7] "MLH"                   "PD-L1(28-8)"           "PD-L1(22C3)"          
    ## [10] "TMB"                   "MSI"                   "MMRD"                 
    ## [13] "Best.Response.(Cycle)" "sex"                   "age"                  
    ## [16] "ATM_mut"               "Responder"             "DCB"

``` r
summarise(df_surv, 
          n = n(), 
          n_PD = sum(PD), 
          n_Death = sum(Death))
```

    ## # A tibble: 2 × 4
    ##   ATM_mut              n  n_PD n_Death
    ##   <chr>            <int> <dbl>   <dbl>
    ## 1 ATM_mut_included     1     0       0
    ## 2 Otherwise           27    16      12

``` r
# Cox proportional hazards model for PFS
coxph(Surv(df_surv$PFS, df_surv$PD) ~ EBV + `PD-L1(28-8)` + `PD-L1(22C3)` + MSI + MMRD + Responder + DCB, data = df_surv)
```

    ## Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    ## Ran out of iterations and did not converge

    ## Call:
    ## coxph(formula = Surv(df_surv$PFS, df_surv$PD) ~ EBV + `PD-L1(28-8)` + 
    ##     `PD-L1(22C3)` + MSI + MMRD + Responder + DCB, data = df_surv)
    ## 
    ##                                               coef
    ## EBVneg                    20.011115346960590244407
    ## EBVNot applicable         20.931316922038792682770
    ## `PD-L1(28-8)`              0.135790199993481353058
    ## `PD-L1(22C3)`              1.350157728138092094738
    ## MSIIndeterminate           0.000000000000000000000
    ## MSIHigh                    0.000000000000000000000
    ## MMRDMMRd                   0.000000000000000000000
    ## Responder                -41.029758050146725167906
    ## DCB                        0.000000000000000000000
    ##                                          exp(coef)
    ## EBVneg             490588057.553904533386230468750
    ## EBVNot applicable 1231276067.793943166732788085937
    ## `PD-L1(28-8)`              1.145441554738545697489
    ## `PD-L1(22C3)`              3.858034003229118269473
    ## MSIIndeterminate           1.000000000000000000000
    ## MSIHigh                    1.000000000000000000000
    ## MMRDMMRd                   1.000000000000000000000
    ## Responder                  0.000000000000000001517
    ## DCB                        1.000000000000000000000
    ##                                           se(coef)      z                   p
    ## EBVneg                     1.400187353659966271735 14.292 <0.0000000000000002
    ## EBVNot applicable      23137.745447469431383069605  0.001               0.999
    ## `PD-L1(28-8)`              0.150362218532861607878  0.903               0.366
    ## `PD-L1(22C3)`              0.093345823595878521517 14.464 <0.0000000000000002
    ## MSIIndeterminate           0.000000000000000000000    NaN                 NaN
    ## MSIHigh                    0.000000000000000000000    NaN                 NaN
    ## MMRDMMRd                   0.000000000000000000000    NaN                 NaN
    ## Responder              16360.856722554175576078705 -0.003               0.998
    ## DCB                    16360.856722554175576078705  0.000               1.000
    ## 
    ## Likelihood ratio test=11.64  on 9 df, p=0.2344
    ## n= 7, number of events= 5 
    ##    (21 observations deleted due to missingness)

``` r
# Cox proportional hazards model for PFS with only PD-L1 related cols.
coxph(Surv(df_surv$PFS, df_surv$PD) ~ `PD-L1(28-8)` + `PD-L1(22C3)`, data = df_surv)
```

    ## Call:
    ## coxph(formula = Surv(df_surv$PFS, df_surv$PD) ~ `PD-L1(28-8)` + 
    ##     `PD-L1(22C3)`, data = df_surv)
    ## 
    ##                   coef exp(coef) se(coef)      z     p
    ## `PD-L1(28-8)`  0.09324   1.09772  0.06295  1.481 0.139
    ## `PD-L1(22C3)` -0.10325   0.90190  0.08057 -1.281 0.200
    ## 
    ## Likelihood ratio test=2.89  on 2 df, p=0.2359
    ## n= 14, number of events= 10 
    ##    (14 observations deleted due to missingness)

``` r
# Cox proportional hazards model for OS
coxph(Surv(df_surv$OS, df_surv$Death) ~ age + sex + HER2 + `PD-L1(28-8)` + `PD-L1(22C3)` + MSI + MMRD + Responder + DCB, data = df_surv)
```

    ## Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    ## Ran out of iterations and did not converge

    ## Call:
    ## coxph(formula = Surv(df_surv$OS, df_surv$Death) ~ age + sex + 
    ##     HER2 + `PD-L1(28-8)` + `PD-L1(22C3)` + MSI + MMRD + Responder + 
    ##     DCB, data = df_surv)
    ## 
    ##                                            coef                      exp(coef)
    ## age                         -5.2762597055427465             0.0051115136488705
    ## sexF                       955.5021924546221044                            Inf
    ## HER2negative                 0.0000000000000000             1.0000000000000000
    ## `PD-L1(28-8)`              -28.0417736884261579             0.0000000000006632
    ## `PD-L1(22C3)`                2.1469479415070620             8.5586968507890706
    ## MSIIndeterminate            27.6717062043878173 1041524065948.1468505859375000
    ## MSIHigh                      0.0000000000000000             1.0000000000000000
    ## MMRDMMRd                     0.0000000000000000             1.0000000000000000
    ## Responder                 -985.2928521589745969             0.0000000000000000
    ## DCB                          0.0000000000000000             1.0000000000000000
    ##                                        se(coef)      z     p
    ## age                         47.1172457589269698 -0.112 0.911
    ## sexF                      2624.4567945064359265  0.364 0.716
    ## HER2negative                 0.0000000000000000    NaN   NaN
    ## `PD-L1(28-8)`               83.9710086272086329 -0.334 0.738
    ## `PD-L1(22C3)`               29.3271639720104673  0.073 0.942
    ## MSIIndeterminate           449.8313250711250930  0.062 0.951
    ## MSIHigh                      0.0000000000000000    NaN   NaN
    ## MMRDMMRd                   449.8313250711250930  0.000 1.000
    ## Responder                 6028.9448483296855557 -0.163 0.870
    ## DCB                       6028.9448483296855557  0.000 1.000
    ## 
    ## Likelihood ratio test=17.05  on 10 df, p=0.07326
    ## n= 10, number of events= 4 
    ##    (18 observations deleted due to missingness)

``` r
# Cox proportional hazards model for OS with only PD-L1 related cols
coxph(Surv(df_surv$OS, df_surv$Death) ~ `PD-L1(28-8)` + `PD-L1(22C3)`, data = df_surv)
```

    ## Call:
    ## coxph(formula = Surv(df_surv$OS, df_surv$Death) ~ `PD-L1(28-8)` + 
    ##     `PD-L1(22C3)`, data = df_surv)
    ## 
    ##                   coef exp(coef) se(coef)      z      p
    ## `PD-L1(28-8)`  0.13080   1.13974  0.07727  1.693 0.0905
    ## `PD-L1(22C3)` -0.03486   0.96574  0.07964 -0.438 0.6616
    ## 
    ## Likelihood ratio test=3.66  on 2 df, p=0.1606
    ## n= 14, number of events= 6 
    ##    (14 observations deleted due to missingness)

# Fit to ATM+TP53+ARID1A group

## Preprocess

``` r
# Create a scatterplot of PD-L1(22C3) vs PD-L1(28-8) with different colors for ATM_mut groups
ggplot(df_surv, aes(x = `PD-L1(22C3)`, y = `PD-L1(28-8)`, color = ATM_mut)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  scale_color_manual(values = c("ATM_mut_included" = "blue", "Otherwise" = "red")) +
  labs(title = "Relationship between PD-L1(22C3) and PD-L1(28-8)",
       subtitle = "Blue points indicate ATM_mut_included samples",
       x = "PD-L1(22C3)",
       y = "PD-L1(28-8)",
       color = "Group") +
  theme_minimal() +
  # Add correlation coefficient
  annotate("text", x = max(df_surv$`PD-L1(22C3)`, na.rm = TRUE) * 0.8, 
           y = max(df_surv$`PD-L1(28-8)`, na.rm = TRUE) * 0.9,
           label = paste("Correlation =", round(cor(df_surv$`PD-L1(22C3)`, df_surv$`PD-L1(28-8)`, 
                                                  use = "complete.obs"), 3)),
           size = 4)
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 14 rows containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning: Removed 14 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-77-1.png)<!-- -->

``` r
which(is.na(df_surv$`PD-L1(28-8)`))
```

    ## [1]  9 10 13 14 23 24 25

``` r
which(is.na(df_surv$`PD-L1(22C3)`))
```

    ## [1]  2  5  7  9 12 19 25 26 28

``` r
# Print indices with NA values in at least one PD-L1 column
both_na_indices <- which(is.na(df_surv$`PD-L1(28-8)`) & is.na(df_surv$`PD-L1(22C3)`))
both_na_indices
```

    ## [1]  9 25

``` r
# Make index column for df_surv
df_surv$index <- 1:nrow(df_surv)
# Print the indices of df_surv that correspond to the NA values in both PD-L1 columns
df_surv$Death[df_surv$index %in% both_na_indices]
```

    ## [1] 0 1

``` r
table(df_surv$Death)
```

    ## 
    ##  0  1 
    ## 16 12

``` r
which(!is.na(df_surv$`PD-L1(28-8)`) & !is.na(df_surv$`PD-L1(22C3)`))
```

    ##  [1]  1  3  4  6  8 11 15 16 17 18 20 21 22 27

``` r
length(which(!is.na(df_surv$`PD-L1(28-8)`) & !is.na(df_surv$`PD-L1(22C3)`)))
```

    ## [1] 14

``` r
# Count rows with NAs in at least one of the two PD-L1 columns
length(which(is.na(df_surv$`PD-L1(28-8)`) | is.na(df_surv$`PD-L1(22C3)`)))
```

    ## [1] 14

## Variable control

By the domain knowledge from Dr. An, 22C3 will be utilized primarily
over 28-8.

``` r
df_surv <- df_surv |> mutate(PDL1 = ifelse(is.na(`PD-L1(22C3)`) == TRUE,
                                                            `PD-L1(28-8)`, `PD-L1(22C3)`))
head(df_surv[, c('PDL1', 'PD-L1(28-8)', 'PD-L1(22C3)')])
```

    ## # A tibble: 6 × 3
    ##    PDL1 `PD-L1(28-8)` `PD-L1(22C3)`
    ##   <dbl>         <dbl>         <dbl>
    ## 1    20             5            20
    ## 2     1             1            NA
    ## 3     5             3             5
    ## 4     5            35             5
    ## 5     1             1            NA
    ## 6    20            10            20

``` r
coxph(Surv(df_surv$OS, df_surv$Death) ~ EBV + PDL1 + MSI + MMRD + Responder + DCB, data = df_surv)
```

    ## Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    ## Loglik converged before variable 5,8 ; coefficient may be infinite.

    ## Call:
    ## coxph(formula = Surv(df_surv$OS, df_surv$Death) ~ EBV + PDL1 + 
    ##     MSI + MMRD + Responder + DCB, data = df_surv)
    ## 
    ##                                  coef           exp(coef)            se(coef)
    ## EBVneg                1.0561974718479     2.8754163232364     1.0976084832688
    ## EBVNot applicable     0.3774815204815     1.4586064891997 43795.7647427787597
    ## PDL1                  0.0137025327620     1.0137968427336     0.0200697696824
    ## MSIIndeterminate                   NA                  NA     0.0000000000000
    ## MSIHigh             -20.6605364442995     0.0000000010647 49375.9676283549270
    ## MMRDMMRd                           NA                  NA     0.0000000000000
    ## Responder             0.6239058203183     1.8662028786167 43795.7647391296268
    ## DCB                 -21.6582304897809     0.0000000003926 37978.4026149316560
    ##                        z     p
    ## EBVneg             0.962 0.336
    ## EBVNot applicable  0.000 1.000
    ## PDL1               0.683 0.495
    ## MSIIndeterminate      NA    NA
    ## MSIHigh            0.000 1.000
    ## MMRDMMRd              NA    NA
    ## Responder          0.000 1.000
    ## DCB               -0.001 1.000
    ## 
    ## Likelihood ratio test=11.57  on 6 df, p=0.0723
    ## n= 16, number of events= 7 
    ##    (12 observations deleted due to missingness)

``` r
table(df_surv$EBV)
```

    ## 
    ##            pos            neg Not applicable 
    ##              8             12              1

``` r
table(df_surv$MSI)
```

    ## 
    ##           MSS Indeterminate          High 
    ##            18             1             4

- ‘Not applicable’ in `EBV` means that the experiment was not conducted.
- ‘MSS’ in `MSI` is the abbreviation of Macrosatellite Stable, so this
  means the corresponding sample shows low level of MSI(microsatellite
  instability) in the test.
- ‘Indeterminate’ in `MSI` means that the corresponding sample shows
  middle level of MSI.

We will use `MMRD` instead of `MSI`.

``` r
df_surv_ordered <- df_surv |> 
  select(-c('PD-L1(28-8)', 'PD-L1(22C3)', 'MSI')) |>
  mutate(EBV = ifelse(EBV == 'neg', -1,
                      ifelse(EBV == 'pos', 1, -1)), 
         ATM_mut = ifelse(ATM_mut == 'ATM_mut_included', 1, 0)) |>  # One-hot-encode the target too.
  mutate(EBV = replace_na(EBV, -1)) # Treat NA and 'Not Applicable' as Negative, said by Dr.An.
table(df_surv_ordered$EBV, useNA = 'ifany')
```

    ## 
    ## -1  1 
    ## 20  8

``` r
skim(df_surv_ordered)
```

|                                                  |                 |
|:-------------------------------------------------|:----------------|
| Name                                             | df_surv_ordered |
| Number of rows                                   | 28              |
| Number of columns                                | 17              |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |                 |
| Column type frequency:                           |                 |
| character                                        | 2               |
| factor                                           | 4               |
| numeric                                          | 10              |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |                 |
| Group variables                                  | ATM_mut         |

Data summary

**Variable type: character**

| skim_variable | ATM_mut | n_missing | complete_rate | min | max | empty | n_unique | whitespace |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| MLH | 0 | 1 | 0.96 | 6 | 45 | 0 | 2 | 0 |
| MLH | 1 | 0 | 1.00 | 6 | 6 | 0 | 1 | 0 |
| Best.Response.(Cycle) | 0 | 0 | 1.00 | 2 | 7 | 0 | 6 | 0 |
| Best.Response.(Cycle) | 1 | 0 | 1.00 | 7 | 7 | 0 | 1 | 0 |

**Variable type: factor**

| skim_variable | ATM_mut | n_missing | complete_rate | ordered | n_unique | top_counts      |
|:--------------|--------:|----------:|--------------:|:--------|---------:|:----------------|
| HER2          |       0 |         0 |          1.00 | FALSE   |        2 | neu: 26, neg: 1 |
| HER2          |       1 |         0 |          1.00 | FALSE   |        1 | neu: 1, neg: 0  |
| TMB           |       0 |         7 |          0.74 | FALSE   |        2 | Low: 12, Hig: 8 |
| TMB           |       1 |         0 |          1.00 | FALSE   |        1 | Low: 1, Hig: 0  |
| MMRD          |       0 |         0 |          1.00 | FALSE   |        2 | MMR: 21, MMR: 6 |
| MMRD          |       1 |         0 |          1.00 | FALSE   |        1 | MMR: 1, MMR: 0  |
| sex           |       0 |         0 |          1.00 | FALSE   |        2 | M: 21, F: 6     |
| sex           |       1 |         0 |          1.00 | FALSE   |        1 | F: 1, M: 0      |

**Variable type: numeric**

| skim_variable | ATM_mut | n_missing | complete_rate | mean | sd | p0 | p25 | p50 | p75 | p100 | hist |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| PFS | 0 | 0 | 1.00 | 4.78 | 3.80 | 1.15 | 1.86 | 3.65 | 5.57 | 13.60 | ▇▅▁▁▂ |
| PFS | 1 | 0 | 1.00 | 11.96 | NA | 11.96 | 11.96 | 11.96 | 11.96 | 11.96 | ▁▁▇▁▁ |
| PD | 0 | 0 | 1.00 | 0.59 | 0.50 | 0.00 | 0.00 | 1.00 | 1.00 | 1.00 | ▆▁▁▁▇ |
| PD | 1 | 0 | 1.00 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| OS | 0 | 0 | 1.00 | 6.37 | 3.27 | 1.61 | 4.40 | 5.78 | 7.80 | 13.60 | ▅▇▆▁▂ |
| OS | 1 | 0 | 1.00 | 11.96 | NA | 11.96 | 11.96 | 11.96 | 11.96 | 11.96 | ▁▁▇▁▁ |
| Death | 0 | 0 | 1.00 | 0.44 | 0.51 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▆ |
| Death | 1 | 0 | 1.00 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| EBV | 0 | 0 | 1.00 | -0.41 | 0.93 | -1.00 | -1.00 | -1.00 | 1.00 | 1.00 | ▇▁▁▁▃ |
| EBV | 1 | 0 | 1.00 | -1.00 | NA | -1.00 | -1.00 | -1.00 | -1.00 | -1.00 | ▁▁▇▁▁ |
| age | 0 | 0 | 1.00 | 62.67 | 11.59 | 25.00 | 56.00 | 64.00 | 72.50 | 77.00 | ▁▁▆▇▇ |
| age | 1 | 0 | 1.00 | 58.00 | NA | 58.00 | 58.00 | 58.00 | 58.00 | 58.00 | ▁▁▇▁▁ |
| Responder | 0 | 0 | 1.00 | 0.30 | 0.47 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▃ |
| Responder | 1 | 0 | 1.00 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |
| DCB | 0 | 0 | 1.00 | 0.37 | 0.49 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▅ |
| DCB | 1 | 0 | 1.00 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |
| index | 0 | 0 | 1.00 | 14.56 | 8.38 | 1.00 | 7.50 | 15.00 | 21.50 | 28.00 | ▇▇▇▇▇ |
| index | 1 | 0 | 1.00 | 13.00 | NA | 13.00 | 13.00 | 13.00 | 13.00 | 13.00 | ▁▁▇▁▁ |
| PDL1 | 0 | 2 | 0.93 | 12.44 | 20.26 | 1.00 | 3.00 | 5.00 | 10.00 | 90.00 | ▇▁▁▁▁ |
| PDL1 | 1 | 0 | 1.00 | 5.00 | NA | 5.00 | 5.00 | 5.00 | 5.00 | 5.00 | ▁▁▇▁▁ |

``` r
# Drop an unimportant factor varible: TMB
df_surv_ordered <- df_surv_ordered |> select(-c('TMB', 'MLH'))
# See the left factor variables
table(df_surv_ordered$HER2)
```

    ## 
    ##  neutral negative 
    ##       27        1

``` r
table(df_surv_ordered$MMRD)
```

    ## 
    ## MMRp MMRd 
    ##   22    6

``` r
table(df_surv_ordered$sex)
```

    ## 
    ##  M  F 
    ## 21  7

In `MMRD` column, - MMRd: mismatch repair deficiency -\> Can lead to
increased mutation rates potentially cancer development - MMRp: mismatch
repair proficiency

``` r
df_surv_ordered <- df_surv_ordered |> 
  mutate(HER2 = ifelse(HER2 == 'neutral', 0, -1),
         MMRD = ifelse(MMRD == 'MMRd', 1, -1),
         sex_fem = ifelse(sex == 'M', 0, 1)) |> 
  select(-c('sex', 'Best.Response.(Cycle)', 'index'))
skim(df_surv_ordered)
```

|                                                  |                 |
|:-------------------------------------------------|:----------------|
| Name                                             | df_surv_ordered |
| Number of rows                                   | 28              |
| Number of columns                                | 13              |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |                 |
| Column type frequency:                           |                 |
| numeric                                          | 12              |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |                 |
| Group variables                                  | ATM_mut         |

Data summary

**Variable type: numeric**

| skim_variable | ATM_mut | n_missing | complete_rate | mean | sd | p0 | p25 | p50 | p75 | p100 | hist |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| PFS | 0 | 0 | 1.00 | 4.78 | 3.80 | 1.15 | 1.86 | 3.65 | 5.57 | 13.60 | ▇▅▁▁▂ |
| PFS | 1 | 0 | 1.00 | 11.96 | NA | 11.96 | 11.96 | 11.96 | 11.96 | 11.96 | ▁▁▇▁▁ |
| PD | 0 | 0 | 1.00 | 0.59 | 0.50 | 0.00 | 0.00 | 1.00 | 1.00 | 1.00 | ▆▁▁▁▇ |
| PD | 1 | 0 | 1.00 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| OS | 0 | 0 | 1.00 | 6.37 | 3.27 | 1.61 | 4.40 | 5.78 | 7.80 | 13.60 | ▅▇▆▁▂ |
| OS | 1 | 0 | 1.00 | 11.96 | NA | 11.96 | 11.96 | 11.96 | 11.96 | 11.96 | ▁▁▇▁▁ |
| Death | 0 | 0 | 1.00 | 0.44 | 0.51 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▆ |
| Death | 1 | 0 | 1.00 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| HER2 | 0 | 0 | 1.00 | -0.04 | 0.19 | -1.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▁▁▇ |
| HER2 | 1 | 0 | 1.00 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| EBV | 0 | 0 | 1.00 | -0.41 | 0.93 | -1.00 | -1.00 | -1.00 | 1.00 | 1.00 | ▇▁▁▁▃ |
| EBV | 1 | 0 | 1.00 | -1.00 | NA | -1.00 | -1.00 | -1.00 | -1.00 | -1.00 | ▁▁▇▁▁ |
| MMRD | 0 | 0 | 1.00 | -0.56 | 0.85 | -1.00 | -1.00 | -1.00 | -1.00 | 1.00 | ▇▁▁▁▂ |
| MMRD | 1 | 0 | 1.00 | -1.00 | NA | -1.00 | -1.00 | -1.00 | -1.00 | -1.00 | ▁▁▇▁▁ |
| age | 0 | 0 | 1.00 | 62.67 | 11.59 | 25.00 | 56.00 | 64.00 | 72.50 | 77.00 | ▁▁▆▇▇ |
| age | 1 | 0 | 1.00 | 58.00 | NA | 58.00 | 58.00 | 58.00 | 58.00 | 58.00 | ▁▁▇▁▁ |
| Responder | 0 | 0 | 1.00 | 0.30 | 0.47 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▃ |
| Responder | 1 | 0 | 1.00 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |
| DCB | 0 | 0 | 1.00 | 0.37 | 0.49 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▅ |
| DCB | 1 | 0 | 1.00 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |
| PDL1 | 0 | 2 | 0.93 | 12.44 | 20.26 | 1.00 | 3.00 | 5.00 | 10.00 | 90.00 | ▇▁▁▁▁ |
| PDL1 | 1 | 0 | 1.00 | 5.00 | NA | 5.00 | 5.00 | 5.00 | 5.00 | 5.00 | ▁▁▇▁▁ |
| sex_fem | 0 | 0 | 1.00 | 0.22 | 0.42 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | ▇▁▁▁▂ |
| sex_fem | 1 | 0 | 1.00 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |

``` r
table(df_surv_ordered$EBV, df_surv_ordered$ATM_mut, useNA = 'ifany')
```

    ##     
    ##       0  1
    ##   -1 19  1
    ##   1   8  0

``` r
# Remove rows that have NA in PDL1, EBV, MSI
nrow(df_surv_ordered)
```

    ## [1] 28

``` r
df_surv_ordered <- df_surv_ordered[!is.na(df_surv_ordered$PDL1), ]
df_surv_ordered <- df_surv_ordered[!is.na(df_surv_ordered$EBV), ]
nrow(df_surv_ordered)
```

    ## [1] 26

``` r
skim(df_surv_ordered)
```

|                                                  |                 |
|:-------------------------------------------------|:----------------|
| Name                                             | df_surv_ordered |
| Number of rows                                   | 26              |
| Number of columns                                | 13              |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_   |                 |
| Column type frequency:                           |                 |
| numeric                                          | 12              |
| \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_ |                 |
| Group variables                                  | ATM_mut         |

Data summary

**Variable type: numeric**

| skim_variable | ATM_mut | n_missing | complete_rate | mean | sd | p0 | p25 | p50 | p75 | p100 | hist |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|
| PFS | 0 | 0 | 1 | 4.60 | 3.54 | 1.15 | 1.87 | 3.65 | 5.55 | 13.60 | ▇▅▁▁▁ |
| PFS | 1 | 0 | 1 | 11.96 | NA | 11.96 | 11.96 | 11.96 | 11.96 | 11.96 | ▁▁▇▁▁ |
| PD | 0 | 0 | 1 | 0.60 | 0.50 | 0.00 | 0.00 | 1.00 | 1.00 | 1.00 | ▅▁▁▁▇ |
| PD | 1 | 0 | 1 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| OS | 0 | 0 | 1 | 6.31 | 2.99 | 1.94 | 4.80 | 5.78 | 7.79 | 13.60 | ▅▇▆▁▂ |
| OS | 1 | 0 | 1 | 11.96 | NA | 11.96 | 11.96 | 11.96 | 11.96 | 11.96 | ▁▁▇▁▁ |
| Death | 0 | 0 | 1 | 0.44 | 0.51 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▆ |
| Death | 1 | 0 | 1 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| HER2 | 0 | 0 | 1 | -0.04 | 0.20 | -1.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▁▁▇ |
| HER2 | 1 | 0 | 1 | 0.00 | NA | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | ▁▁▇▁▁ |
| EBV | 0 | 0 | 1 | -0.44 | 0.92 | -1.00 | -1.00 | -1.00 | 1.00 | 1.00 | ▇▁▁▁▃ |
| EBV | 1 | 0 | 1 | -1.00 | NA | -1.00 | -1.00 | -1.00 | -1.00 | -1.00 | ▁▁▇▁▁ |
| MMRD | 0 | 0 | 1 | -0.60 | 0.82 | -1.00 | -1.00 | -1.00 | -1.00 | 1.00 | ▇▁▁▁▂ |
| MMRD | 1 | 0 | 1 | -1.00 | NA | -1.00 | -1.00 | -1.00 | -1.00 | -1.00 | ▁▁▇▁▁ |
| age | 0 | 0 | 1 | 62.92 | 12.02 | 25.00 | 56.00 | 64.00 | 73.00 | 77.00 | ▁▁▆▆▇ |
| age | 1 | 0 | 1 | 58.00 | NA | 58.00 | 58.00 | 58.00 | 58.00 | 58.00 | ▁▁▇▁▁ |
| Responder | 0 | 0 | 1 | 0.32 | 0.48 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▃ |
| Responder | 1 | 0 | 1 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |
| DCB | 0 | 0 | 1 | 0.36 | 0.49 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | ▇▁▁▁▅ |
| DCB | 1 | 0 | 1 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |
| PDL1 | 0 | 0 | 1 | 12.44 | 20.26 | 1.00 | 3.00 | 5.00 | 10.00 | 90.00 | ▇▁▁▁▁ |
| PDL1 | 1 | 0 | 1 | 5.00 | NA | 5.00 | 5.00 | 5.00 | 5.00 | 5.00 | ▁▁▇▁▁ |
| sex_fem | 0 | 0 | 1 | 0.24 | 0.44 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | ▇▁▁▁▂ |
| sex_fem | 1 | 0 | 1 | 1.00 | NA | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | ▁▁▇▁▁ |

### Cox with combined PD-L1

This Shows that MMRD and ATM_mut are beneficial predictors for PFS, and
age is a beneficial predictors for OS.

``` r
coxph(Surv(df_surv_ordered$PFS, df_surv_ordered$PD) ~ EBV + PDL1 + MMRD + ATM_mut + HER2 + age + sex_fem, data = df_surv_ordered)
```

    ## Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    ## Ran out of iterations and did not converge

    ## Call:
    ## coxph(formula = Surv(df_surv_ordered$PFS, df_surv_ordered$PD) ~ 
    ##     EBV + PDL1 + MMRD + ATM_mut + HER2 + age + sex_fem, data = df_surv_ordered)
    ## 
    ##                     coef        exp(coef)         se(coef)      z      p
    ## EBV       -0.68663376049    0.50326733926    0.33583592175 -2.045 0.0409
    ## PDL1      -0.00895680720    0.99108318550    0.01565754269 -0.572 0.5673
    ## MMRD      -0.99302764465    0.37045339018    0.52138696960 -1.905 0.0568
    ## ATM_mut  -17.94910857993    0.00000001603 5532.61071208307 -0.003 0.9974
    ## HER2      -1.90039816451    0.14950907816    1.41435440103 -1.344 0.1791
    ## age       -0.02405171247    0.97623522492    0.02188699755 -1.099 0.2718
    ## sex_fem   -0.34225033607    0.71017040119    0.62222297253 -0.550 0.5823
    ## 
    ## Likelihood ratio test=14.93  on 7 df, p=0.03689
    ## n= 26, number of events= 15

``` r
coxph(Surv(df_surv_ordered$OS, df_surv_ordered$Death) ~ EBV + PDL1 + MMRD + ATM_mut + HER2 + age + sex_fem, data = df_surv_ordered)
```

    ## Warning in coxph.fit(X, Y, istrat, offset, init, control, weights = weights, :
    ## Loglik converged before variable 4 ; coefficient may be infinite.

    ## Call:
    ## coxph(formula = Surv(df_surv_ordered$OS, df_surv_ordered$Death) ~ 
    ##     EBV + PDL1 + MMRD + ATM_mut + HER2 + age + sex_fem, data = df_surv_ordered)
    ## 
    ##                       coef          exp(coef)           se(coef)      z     p
    ## EBV        -0.576630341228     0.561788215292     0.490973911567 -1.174 0.240
    ## PDL1        0.003956687429     1.003964525451     0.017283395532  0.229 0.819
    ## MMRD       -0.188271913745     0.828389426134     0.569142119393 -0.331 0.741
    ## ATM_mut   -18.922021414914     0.000000006057 11570.945008313585 -0.002 0.999
    ## HER2                    NA                 NA     0.000000000000     NA    NA
    ## age        -0.049539007641     0.951668035088     0.038000364416 -1.304 0.192
    ## sex_fem    -0.063261247686     0.938698208918     0.917041372682 -0.069 0.945
    ## 
    ## Likelihood ratio test=12.17  on 6 df, p=0.05823
    ## n= 26, number of events= 11

## Model Fitting

### Target: Responder

#### Logit.Reg.+LASSO

``` r
# Apply logistic regression with LASSO for 'Responder' prediction
library(glmnet)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## Loaded glmnet 4.1-9

``` r
# Prepare predictor matrix and response vector
x <- as.matrix(df_surv_ordered[, !names(df_surv_ordered) %in% c("Responder", "PD", "PFS", "Death", "OS", "HER2",
                                                                "DCB")])    # Might be weird to use post-knowledge
y <- as.numeric(df_surv_ordered$Responder)

# alpha = 1: LASSO, 0.5: Elastic Net
fit <- glmnet(x, y, family = "gaussian", # gaussian means logistic regression on regression problem
              type.measure = "mse",
              standardize = TRUE, # standardize the data -> Can differ the coefficient's axis
              # standardize = FALSE,
              alpha = 1)  # alpha=1 means LASSO, 0.5 means Elastic Net

# Extract coefficients by glmnet::predict
beta_pred <- predict(fit, type = "coefficients")
lambda_seq <- fit$lambda

# Make coef.s tidy
coef_list <- lapply(1:length(lambda_seq), function(i) {
  beta <- as.matrix(beta_pred[, i])  # ith coefficient for lambda
  df <- data.frame(variable = rownames(beta),
                   coefficient = as.numeric(beta),
                   lambda = lambda_seq[i])
  df
})

coef_long <- bind_rows(coef_list) |>
  filter(variable != "(Intercept)")  # Exclude the intercept

ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(linewidth = 1.2) +  # Increased line thickness
  theme_minimal() +
  labs(
      # title = "Elastic Net Regularization Path",
      title = "LASSO Regularization Path of predicting 'Responder' in Logit.Reg.",
       x = "log(lambda)",
       y = "Coefficient") +
  theme(legend.position = "right")
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-98-1.png)<!-- -->

#### Only logit (No LASSO)

``` r
# Fit logistic reg. with scaled x
glm_fit <- glm(y ~ scale(x), family = binomial())

coef_summary <- summary(glm_fit)$coefficients
rownames(coef_summary) <- gsub("scale(x)", "", rownames(coef_summary), fixed = TRUE)
coef_summary
```

    ##               Estimate   Std. Error      z value  Pr(>|z|)
    ## (Intercept) -4.7538426 1021.1320535 -0.004655463 0.9962855
    ## EBV          9.2539388 1568.5961782  0.005899504 0.9952929
    ## MMRD         7.6713864 1393.7342526  0.005504196 0.9956083
    ## age          0.8809583    0.8873556  0.992790674 0.3208120
    ## ATM_mut      4.1210365 2109.0354752  0.001953991 0.9984409
    ## PDL1         0.3029103    0.5659147  0.535257811 0.5924716
    ## sex_fem      8.5141505 1568.5962207  0.005427879 0.9956692

``` r
coef_df <- as.data.frame(coef_summary)
coef_df$variable <- rownames(coef_df)
# Store rows except of 'intercept'
coef_df <- coef_df[coef_df$variable != "(Intercept)", ]

# Sort variables by descending order of their absolutes
coef_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  arrange(desc(abs_coef))

ggplot(coef_df, aes(x = reorder(variable, abs_coef), y = Estimate)) +
  geom_col(aes(fill = Estimate > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance (Logistic Regression Coefficients)",
       x = "Variable",
       y = "Coefficient") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Direction"))
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-100-1.png)<!-- -->
\#### Random Forest

``` r
library(randomForest)
```

    ## randomForest 4.7-1.2

    ## Type rfNews() to see new features/changes/bug fixes.

    ## 
    ## Attaching package: 'randomForest'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     margin

``` r
# Fit a Random Forest model
rf <- randomForest(y = factor(y), x = scale(x), importance = TRUE)

# Display variable importance
importance(rf)
```

    ##                  0          1 MeanDecreaseAccuracy MeanDecreaseGini
    ## EBV      9.0220056 10.6712577           11.7558244        2.0920587
    ## MMRD    -5.0625948  0.5418028           -4.2043250        0.4899436
    ## age     -0.4303729  1.8035404            0.7863433        3.4305492
    ## ATM_mut  0.0000000  0.0000000            0.0000000        0.6113319
    ## PDL1     2.1195006  1.3477472            2.4126160        2.2067797
    ## sex_fem -3.0106767 -0.3103555           -2.8392402        0.4893842

``` r
varImpPlot(rf)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-102-1.png)<!-- -->
\#### Combine LASSO with Logit&RF

``` r
# Extract importance measures
rf_importance <- as.data.frame(importance(rf))
rf_importance$variable <- rownames(rf_importance)

# Merge all information by variable
merged_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  left_join(rf_importance, by = "variable")
merged_df
```

    ##    Estimate   Std. Error     z value  Pr(>|z|) variable  abs_coef          0
    ## 1 9.2539388 1568.5961782 0.005899504 0.9952929      EBV 9.2539388  9.0220056
    ## 2 8.5141505 1568.5962207 0.005427879 0.9956692  sex_fem 8.5141505 -3.0106767
    ## 3 7.6713864 1393.7342526 0.005504196 0.9956083     MMRD 7.6713864 -5.0625948
    ## 4 4.1210365 2109.0354752 0.001953991 0.9984409  ATM_mut 4.1210365  0.0000000
    ## 5 0.8809583    0.8873556 0.992790674 0.3208120      age 0.8809583 -0.4303729
    ## 6 0.3029103    0.5659147 0.535257811 0.5924716     PDL1 0.3029103  2.1195006
    ##            1 MeanDecreaseAccuracy MeanDecreaseGini
    ## 1 10.6712577           11.7558244        2.0920587
    ## 2 -0.3103555           -2.8392402        0.4893842
    ## 3  0.5418028           -4.2043250        0.4899436
    ## 4  0.0000000            0.0000000        0.6113319
    ## 5  1.8035404            0.7863433        3.4305492
    ## 6  1.3477472            2.4126160        2.2067797

``` r
# ---- SELECT SORTING KEY ----
# Default
sorting_key <- merged_df$abs_coef
# sorting_key <- merged_df$MeanDecreaseAccuracy
# sorting_key <- merged_df$MeanDecreaseGini

# ---- Sort ----
merged_df$variable <- factor(merged_df$variable, levels = merged_df$variable[order(sorting_key, decreasing = FALSE)])

# ---- Adjust legend's location in Logit.Reg. ----
legend_pos <- c(0.7, 0.2)  # ex) "bottom", "left", c(0.9, 0.1)

# ---- Plot: Logistic Regression ----
p_logit <- ggplot(merged_df, aes(x = variable, y = Estimate)) +
  geom_col(aes(fill = Estimate > 0), width = 0.6) +
  coord_flip() +
  # geom_text(aes(label = sprintf("%.2f", Estimate)), 
  #           vjust = ifelse(merged_df$Estimate > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Positively\ncorrelated to\nResponder == 1")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = legend_pos,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Coefficient",
       y = "Estimate")
```

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## ℹ Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# ---- Plot: RF Accuracy (Modified to resemble varImpPlot) ----
p_mda <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseAccuracy), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Accuracy",
       y = "MeanDecreaseAccuracy")

# ---- Plot: RF Gini (Modified to resemble varImpPlot) ----
p_gini <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseGini)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseGini), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Gini",
       y = "MeanDecreaseGini")


# ---- Combine all ----
library(grid)
library(gridExtra)

# Text labels
title_all <- textGrob("Scaled Variables' indices related to Target", gp = gpar(fontsize = 14, fontface = "bold"),
                      just = 'center') # center-align in entire plot
predictor_var_text <- textGrob("Target: Responder(=1; Non-responder=0)", gp = gpar(fontsize = 12),
                               just = 'center')
title_logit <- textGrob("Logistic Regression", gp = gpar(fontsize = 13, fontface = "bold"))
title_rf    <- textGrob("Random Forest", gp = gpar(fontsize = 13, fontface = "bold"))



# subtitle row: 1st column = title_logit, 2nd+3rd column merged = title_rf
subtitle_row <- arrangeGrob(
  grobs = list(title_logit, title_rf),
  ncol = 3,
  layout_matrix = rbind(c(1, 2, 2))
)

# Plot row: actual plots
plot_row <- arrangeGrob(
  grobs = list(p_logit, p_mda, p_gini),
  ncol = 3,
  widths = c(2.5, 2.5, 2.5)
)

# Combine titles and plots
grid.arrange(
  title_all,
  predictor_var_text,
  subtitle_row,
  plot_row,
  nrow = 4,
  heights = c(0.12, 0.08, 0.12, 1)  # control height of title vs plot area
)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-104-1.png)<!-- -->
\### Target: ATM_mut \#### Logit.Reg.+LASSO

``` r
# Apply logistic regression with LASSO for 'ATM_mut' prediction
library(glmnet)

# Prepare predictor matrix and response vector
x <- as.matrix(df_surv_ordered[, !names(df_surv_ordered) %in% c("ATM_mut", "PD", "PFS", "Death", "OS", "HER2",
                                                                "Responder")])    # Might be weird to use post-knowledge
y <- as.numeric(df_surv_ordered$ATM_mut)

# alpha = 1: LASSO, 0.5: Elastic Net
fit <- glmnet(x, y, family = "gaussian", # gaussian means logistic regression on regression problem
              type.measure = "mse",
              standardize = TRUE, # standardize the data -> Can differ the coefficient's axis
              # standardize = FALSE,
              alpha = 1)  # alpha=1 means LASSO, 0.5 means Elastic Net

# 계수 추출 (glmnet predict를 이용)
beta_pred <- predict(fit, type = "coefficients")
lambda_seq <- fit$lambda

# 계수 tidy하게 변환
coef_list <- lapply(1:length(lambda_seq), function(i) {
  beta <- as.matrix(beta_pred[, i])  # i번째 lambda에 대한 계수
  df <- data.frame(variable = rownames(beta),
                   coefficient = as.numeric(beta),
                   lambda = lambda_seq[i])
  df
})

coef_long <- bind_rows(coef_list) |>
  filter(variable != "(Intercept)")  # Exclude the intercept

ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(linewidth = 1.2) +  # Increased line thickness
  theme_minimal() +
  labs(
      # title = "Elastic Net Regularization Path",
      title = "LASSO Regularization Path of predicting ATM_mut in Logit.Reg.",
       x = "log(lambda)",
       y = "Coefficient") +
  theme(legend.position = "right")
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-105-1.png)<!-- -->
\#### Only logit (No LASSO)

``` r
# Fit logistic reg. with scaled x
glm_fit <- glm(y ~ scale(x), family = binomial())
```

    ## Warning: glm.fit: algorithm did not converge

    ## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred

``` r
coef_summary <- summary(glm_fit)$coefficients
rownames(coef_summary) <- gsub("scale(x)", "", rownames(coef_summary), fixed = TRUE)
coef_summary
```

    ##               Estimate Std. Error       z value  Pr(>|z|)
    ## (Intercept) -128.25848  100047.68 -0.0012819736 0.9989771
    ## EBV          -55.75433   57750.54 -0.0009654340 0.9992297
    ## MMRD         -11.84778   57222.69 -0.0002070468 0.9998348
    ## age           50.81602   47209.66  0.0010763903 0.9991412
    ## DCB           42.21893   38781.79  0.0010886276 0.9991314
    ## PDL1          20.07951   73321.22  0.0002738567 0.9997815
    ## sex_fem       53.24791   48469.02  0.0010985966 0.9991234

``` r
coef_df <- as.data.frame(coef_summary)
coef_df$variable <- rownames(coef_df)
# Store rows except of 'intercept'
coef_df <- coef_df[coef_df$variable != "(Intercept)", ]

# Sort variables by descending order of their absolutes
coef_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  arrange(desc(abs_coef))

ggplot(coef_df, aes(x = reorder(variable, abs_coef), y = Estimate)) +
  geom_col(aes(fill = Estimate > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance (Logistic Regression Coefficients)",
       x = "Variable",
       y = "Coefficient") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Direction"))
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-107-1.png)<!-- -->
\#### Random Forest

``` r
library(randomForest)
# Fit a Random Forest model
rf <- randomForest(y = factor(y), x = scale(x), importance = TRUE)

# Display variable importance
importance(rf)
```

    ##                  0 1 MeanDecreaseAccuracy MeanDecreaseGini
    ## EBV      0.6128955 0            0.6128955       0.06709730
    ## MMRD     0.1530560 0            0.1530560       0.01784055
    ## age      0.1133665 0            0.1133665       1.05348381
    ## DCB      4.2324398 0            4.2324398       0.47734480
    ## PDL1    -2.4581029 0           -2.4581029       0.16864535
    ## sex_fem -3.5439419 0           -3.5439419       0.39849882

``` r
varImpPlot(rf)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-109-1.png)<!-- -->
\#### Combine LASSO with Logit&RF

``` r
# Extract importance measures
rf_importance <- as.data.frame(importance(rf))
rf_importance$variable <- rownames(rf_importance)

# Merge all information by variable
merged_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  left_join(rf_importance, by = "variable")
merged_df
```

    ##    Estimate Std. Error       z value  Pr(>|z|) variable abs_coef          0 1
    ## 1 -55.75433   57750.54 -0.0009654340 0.9992297      EBV 55.75433  0.6128955 0
    ## 2  53.24791   48469.02  0.0010985966 0.9991234  sex_fem 53.24791 -3.5439419 0
    ## 3  50.81602   47209.66  0.0010763903 0.9991412      age 50.81602  0.1133665 0
    ## 4  42.21893   38781.79  0.0010886276 0.9991314      DCB 42.21893  4.2324398 0
    ## 5  20.07951   73321.22  0.0002738567 0.9997815     PDL1 20.07951 -2.4581029 0
    ## 6 -11.84778   57222.69 -0.0002070468 0.9998348     MMRD 11.84778  0.1530560 0
    ##   MeanDecreaseAccuracy MeanDecreaseGini
    ## 1            0.6128955       0.06709730
    ## 2           -3.5439419       0.39849882
    ## 3            0.1133665       1.05348381
    ## 4            4.2324398       0.47734480
    ## 5           -2.4581029       0.16864535
    ## 6            0.1530560       0.01784055

``` r
# ---- SELECT SORTING KEY ----
# Default
sorting_key <- merged_df$abs_coef
# sorting_key <- merged_df$MeanDecreaseAccuracy
# sorting_key <- merged_df$MeanDecreaseGini

# ---- Sort ----
merged_df$variable <- factor(merged_df$variable, levels = merged_df$variable[order(sorting_key, decreasing = FALSE)])

# ---- Adjust legend's location in Logit.Reg. ----
legend_pos <- c(1, 0.2)  # ex) "bottom", "left", c(0.9, 0.1)

# ---- Plot: Logistic Regression ----
p_logit <- ggplot(merged_df, aes(x = variable, y = Estimate)) +
  geom_col(aes(fill = Estimate > 0), width = 0.6) +
  coord_flip() +
  # geom_text(aes(label = sprintf("%.2f", Estimate)), 
  #           vjust = ifelse(merged_df$Estimate > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Positively\ncorrelated to\nATMmut_inc. == 1")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = legend_pos,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Coefficient",
       y = "Estimate")

# ---- Plot: RF Accuracy (Modified to resemble varImpPlot) ----
p_mda <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseAccuracy), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Accuracy",
       y = "MeanDecreaseAccuracy")

# ---- Plot: RF Gini (Modified to resemble varImpPlot) ----
p_gini <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseGini)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseGini), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Gini",
       y = "MeanDecreaseGini")


# ---- Combine all ----
library(grid)
library(gridExtra)

# Text labels
title_all <- textGrob("Scaled Variables' indices related to Target", gp = gpar(fontsize = 14, fontface = "bold"),
                      just = 'center') # center-align in entire plot
predictor_var_text <- textGrob("Target: ATMmut_included(=1; Otherwise=0)", gp = gpar(fontsize = 12),
                               just = 'center')
title_logit <- textGrob("Logistic Regression", gp = gpar(fontsize = 13, fontface = "bold"))
title_rf    <- textGrob("Random Forest", gp = gpar(fontsize = 13, fontface = "bold"))



# subtitle row: 1st column = title_logit, 2nd+3rd column merged = title_rf
subtitle_row <- arrangeGrob(
  grobs = list(title_logit, title_rf),
  ncol = 3,
  layout_matrix = rbind(c(1, 2, 2))
)

# Plot row: actual plots
plot_row <- arrangeGrob(
  grobs = list(p_logit, p_mda, p_gini),
  ncol = 3,
  widths = c(2.5, 2.5, 2.5)
)

# Combine titles and plots
grid.arrange(
  title_all,
  predictor_var_text,
  subtitle_row,
  plot_row,
  nrow = 4,
  heights = c(0.12, 0.08, 0.12, 1)  # control height of title vs plot area
)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-111-1.png)<!-- -->
\### Target: Alive

``` r
# Create a column for 'Alive'
df_surv_alive_included <- df_surv_ordered |> mutate(Alive = ifelse(Death == 1, 0, 1))
table(df_surv_alive_included$Alive)
```

    ## 
    ##  0  1 
    ## 11 15

``` r
table(df_surv_ordered$Death)
```

    ## 
    ##  0  1 
    ## 15 11

#### Logit.Reg.+LASSO

``` r
# Apply logistic regression with LASSO for 'Death' prediction
library(glmnet)

# Prepare predictor matrix and response vector
x <- as.matrix(df_surv_alive_included[, !names(df_surv_alive_included) %in% c("Responder", "PD", "PFS", "Alive", "Death", "OS", "HER2",
                                                                "DCB")])    # Might be weird to use post-knowledge
y <- as.numeric(df_surv_alive_included$Alive)

# alpha = 1: LASSO, 0.5: Elastic Net
fit <- glmnet(x, y, family = "gaussian", # gaussian means logistic regression on regression problem
              type.measure = "mse",
              standardize = TRUE, # standardize the data -> Can differ the coefficient's axis
              # standardize = FALSE,
              alpha = 1)  # alpha=1 means LASSO, 0.5 means Elastic Net

# 계수 추출 (glmnet predict를 이용)
beta_pred <- predict(fit, type = "coefficients")
lambda_seq <- fit$lambda

# 계수 tidy하게 변환
coef_list <- lapply(1:length(lambda_seq), function(i) {
  beta <- as.matrix(beta_pred[, i])  # i번째 lambda에 대한 계수
  df <- data.frame(variable = rownames(beta),
                   coefficient = as.numeric(beta),
                   lambda = lambda_seq[i])
  df
})

coef_long <- bind_rows(coef_list) |>
  filter(variable != "(Intercept)")  # Exclude the intercept

ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(linewidth = 1.2) +  # Increased line thickness
  theme_minimal() +
  labs(
      # title = "Elastic Net Regularization Path",
      title = "LASSO Regularization Path to target: Alive(=1, Death = 0)",
       x = "log(lambda)",
       y = "Coefficient") +
  theme(legend.position = "right")
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-113-1.png)<!-- -->
\#### Only logit (No LASSO)

``` r
# Fit logistic reg. with scaled x
glm_fit <- glm(y ~ scale(x), family = binomial())

coef_summary <- summary(glm_fit)$coefficients
rownames(coef_summary) <- gsub("scale(x)", "", rownames(coef_summary), fixed = TRUE)
coef_summary
```

    ##               Estimate  Std. Error      z value  Pr(>|z|)
    ## (Intercept)  0.8458302  92.2913413  0.009164784 0.9926877
    ## EBV          0.6015625   0.5480986  1.097544325 0.2724035
    ## MMRD         0.6456237   0.5588171  1.155339834 0.2479513
    ## age          1.0785158   0.6713721  1.606435292 0.1081783
    ## ATM_mut      3.4761221 470.5894894  0.007386740 0.9941063
    ## PDL1         0.1210873   0.5338346  0.226825470 0.8205595
    ## sex_fem     -0.1118696   0.5807474 -0.192630367 0.8472485

``` r
coef_df <- as.data.frame(coef_summary)
coef_df$variable <- rownames(coef_df)
# Store rows except of 'intercept'
coef_df <- coef_df[coef_df$variable != "(Intercept)", ]

# Sort variables by descending order of their absolutes
coef_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  arrange(desc(abs_coef))

ggplot(coef_df, aes(x = reorder(variable, abs_coef), y = Estimate)) +
  geom_col(aes(fill = Estimate > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Logistic Regression Coefficients with scaled variables",
       x = "Variable",
       y = "Coefficient") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Direction"))
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-115-1.png)<!-- -->
\#### Random Forest

``` r
library(randomForest)
# Fit a Random Forest model
rf <- randomForest(y = factor(y), x = scale(x), importance = TRUE)

# Display variable importance
importance(rf)
```

    ##                  0          1 MeanDecreaseAccuracy MeanDecreaseGini
    ## EBV     -2.3187811 -3.2526132            -3.884377        0.6682674
    ## MMRD     0.5569237  2.6608725             2.119968        0.8779311
    ## age     -1.7505246 -0.6170811            -1.746829        4.5697328
    ## ATM_mut  0.0000000  0.0000000             0.000000        0.3795344
    ## PDL1    -3.4366617 -4.7728756            -5.332083        1.9811966
    ## sex_fem -5.4508498  0.4779998            -3.488104        0.6562875

``` r
varImpPlot(rf)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-117-1.png)<!-- -->
\#### Combine LASSOwithLogit, RF

``` r
# Extract importance measures
rf_importance <- as.data.frame(importance(rf))
rf_importance$variable <- rownames(rf_importance)

# Merge all information by variable
merged_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  left_join(rf_importance, by = "variable")
merged_df
```

    ##     Estimate  Std. Error     z value  Pr(>|z|) variable  abs_coef          0
    ## 1  3.4761221 470.5894894  0.00738674 0.9941063  ATM_mut 3.4761221  0.0000000
    ## 2  1.0785158   0.6713721  1.60643529 0.1081783      age 1.0785158 -1.7505246
    ## 3  0.6456237   0.5588171  1.15533983 0.2479513     MMRD 0.6456237  0.5569237
    ## 4  0.6015625   0.5480986  1.09754433 0.2724035      EBV 0.6015625 -2.3187811
    ## 5  0.1210873   0.5338346  0.22682547 0.8205595     PDL1 0.1210873 -3.4366617
    ## 6 -0.1118696   0.5807474 -0.19263037 0.8472485  sex_fem 0.1118696 -5.4508498
    ##            1 MeanDecreaseAccuracy MeanDecreaseGini
    ## 1  0.0000000             0.000000        0.3795344
    ## 2 -0.6170811            -1.746829        4.5697328
    ## 3  2.6608725             2.119968        0.8779311
    ## 4 -3.2526132            -3.884377        0.6682674
    ## 5 -4.7728756            -5.332083        1.9811966
    ## 6  0.4779998            -3.488104        0.6562875

``` r
# ---- SELECT SORTING KEY ----
# Default
sorting_key <- merged_df$abs_coef
# sorting_key <- merged_df$MeanDecreaseAccuracy
# sorting_key <- merged_df$MeanDecreaseGini

# ---- Sort ----
merged_df$variable <- factor(merged_df$variable, levels = merged_df$variable[order(sorting_key, decreasing = FALSE)])

# ---- Adjust legend's location in Logit.Reg. ----
legend_pos <- c(0.6, 0.4)  # ex) "bottom", "left", c(0.9, 0.1)

# ---- Plot: Logistic Regression ----
p_logit <- ggplot(merged_df, aes(x = variable, y = Estimate)) +
  geom_col(aes(fill = Estimate > 0), width = 0.6) +
  coord_flip() +
  # geom_text(aes(label = sprintf("%.2f", Estimate)), 
  #           vjust = ifelse(merged_df$Estimate > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Positively\ncorrelated to\nAlive == 1")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = legend_pos,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Coefficient",
       y = "Estimate")

# ---- Plot: RF Accuracy (Modified to resemble varImpPlot) ----
p_mda <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseAccuracy), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Accuracy",
       y = "MeanDecreaseAccuracy")

# ---- Plot: RF Gini (Modified to resemble varImpPlot) ----
p_gini <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseGini)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseGini), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Gini",
       y = "MeanDecreaseGini")


# ---- Combine all ----
library(grid)
library(gridExtra)

# Text labels
title_all <- textGrob("Scaled Variables' indices related to Target", gp = gpar(fontsize = 14, fontface = "bold"),
                      just = 'center') # center-align in entire plot
predictor_var_text <- textGrob("Target: Alive(=1; Death=0)", gp = gpar(fontsize = 12),
                               just = 'center')
title_logit <- textGrob("Logistic Regression", gp = gpar(fontsize = 13, fontface = "bold"))
title_rf    <- textGrob("Random Forest", gp = gpar(fontsize = 13, fontface = "bold"))



# subtitle row: 1st column = title_logit, 2nd+3rd column merged = title_rf
subtitle_row <- arrangeGrob(
  grobs = list(title_logit, title_rf),
  ncol = 3,
  layout_matrix = rbind(c(1, 2, 2))
)

# Plot row: actual plots
plot_row <- arrangeGrob(
  grobs = list(p_logit, p_mda, p_gini),
  ncol = 3,
  widths = c(2.5, 2.5, 2.5)
)

# Combine titles and plots
grid.arrange(
  title_all,
  predictor_var_text,
  subtitle_row,
  plot_row,
  nrow = 4,
  heights = c(0.12, 0.08, 0.12, 1)  # control height of title vs plot area
)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-119-1.png)<!-- -->

### Target: DiseaseDiagnosed

``` r
# Create a column for 'DiseaseDiagnosed'
df_surv_alive_included <- df_surv_ordered |> mutate(DiseaseDiagnosed = ifelse(PD == 1, 0, 1))
table(df_surv_alive_included$DiseaseDiagnosed)
```

    ## 
    ##  0  1 
    ## 15 11

``` r
table(df_surv_ordered$PD)
```

    ## 
    ##  0  1 
    ## 11 15

#### Logit.Reg.+LASSO

``` r
# Apply logistic regression with LASSO for 'Death' prediction
library(glmnet)

# Prepare predictor matrix and response vector
x <- as.matrix(df_surv_alive_included[, !names(df_surv_alive_included) %in% c("Responder", "PD", "PFS", "DiseaseDiagnosed", "Death", "OS", "HER2",
                                                                "DCB")])    # Might be weird to use post-knowledge
y <- as.numeric(df_surv_alive_included$DiseaseDiagnosed)

# alpha = 1: LASSO, 0.5: Elastic Net
fit <- glmnet(x, y, family = "gaussian", # gaussian means logistic regression on regression problem
              type.measure = "mse",
              standardize = TRUE, # standardize the data -> Can differ the coefficient's axis
              # standardize = FALSE,
              alpha = 1)  # alpha=1 means LASSO, 0.5 means Elastic Net

# 계수 추출 (glmnet predict를 이용)
beta_pred <- predict(fit, type = "coefficients")
lambda_seq <- fit$lambda

# 계수 tidy하게 변환
coef_list <- lapply(1:length(lambda_seq), function(i) {
  beta <- as.matrix(beta_pred[, i])  # i번째 lambda에 대한 계수
  df <- data.frame(variable = rownames(beta),
                   coefficient = as.numeric(beta),
                   lambda = lambda_seq[i])
  df
})

coef_long <- bind_rows(coef_list) |>
  filter(variable != "(Intercept)")  # Exclude the intercept

ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(linewidth = 1.2) +  # Increased line thickness
  theme_minimal() +
  labs(
      # title = "Elastic Net Regularization Path",
      title = "LASSO Regularization Path to target: Disease Diagnosed(=1, otherwise = 0)",
       x = "log(lambda)",
       y = "Coefficient") +
  theme(legend.position = "right")
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-121-1.png)<!-- -->
\#### Only logit (No LASSO)

``` r
# Fit logistic reg. with scaled x
glm_fit <- glm(y ~ scale(x), family = binomial())

coef_summary <- summary(glm_fit)$coefficients
rownames(coef_summary) <- gsub("scale(x)", "", rownames(coef_summary), fixed = TRUE)
coef_summary
```

    ##                Estimate  Std. Error       z value   Pr(>|z|)
    ## (Intercept) 0.002735618  92.2917318 0.00002964099 0.99997635
    ## EBV         1.079678527   0.6620937 1.63070348405 0.10295290
    ## MMRD        1.497335717   0.6657342 2.24914950547 0.02450299
    ## age         1.080202324   0.8337743 1.29555719254 0.19512809
    ## ATM_mut     3.531750709 470.5894972 0.00750495013 0.99401197
    ## PDL1        0.381657555   0.5335581 0.71530640573 0.47441973
    ## sex_fem     0.727978233   0.7319215 0.99461241835 0.31992481

``` r
coef_df <- as.data.frame(coef_summary)
coef_df$variable <- rownames(coef_df)
# Store rows except of 'intercept'
coef_df <- coef_df[coef_df$variable != "(Intercept)", ]

# Sort variables by descending order of their absolutes
coef_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  arrange(desc(abs_coef))

ggplot(coef_df, aes(x = reorder(variable, abs_coef), y = Estimate)) +
  geom_col(aes(fill = Estimate > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance (Logistic Regression Coefficients)",
       x = "Variable",
       y = "Coefficient") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Direction"))
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-123-1.png)<!-- -->
\#### Random Forest

``` r
library(randomForest)
# Fit a Random Forest model
rf <- randomForest(y = factor(y), x = scale(x), importance = TRUE)

# Display variable importance
importance(rf)
```

    ##                 0          1 MeanDecreaseAccuracy MeanDecreaseGini
    ## EBV     -1.005123  1.0805885            0.4480967        0.8550596
    ## MMRD     5.686793  8.4643898            8.9514974        1.6613782
    ## age      2.757752  0.4536143            1.7002216        3.5860663
    ## ATM_mut  0.000000  0.0000000            0.0000000        0.5884010
    ## PDL1     2.616437 -0.3073140            1.7497756        2.4945216
    ## sex_fem -4.594646  0.9889172           -2.9122258        0.5571359

``` r
varImpPlot(rf)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-125-1.png)<!-- -->
\#### Combine LASSO with Logit&RF

``` r
# Extract importance measures
rf_importance <- as.data.frame(importance(rf))
rf_importance$variable <- rownames(rf_importance)

# Merge all information by variable
merged_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  left_join(rf_importance, by = "variable")
merged_df
```

    ##    Estimate  Std. Error    z value   Pr(>|z|) variable  abs_coef         0
    ## 1 3.5317507 470.5894972 0.00750495 0.99401197  ATM_mut 3.5317507  0.000000
    ## 2 1.4973357   0.6657342 2.24914951 0.02450299     MMRD 1.4973357  5.686793
    ## 3 1.0802023   0.8337743 1.29555719 0.19512809      age 1.0802023  2.757752
    ## 4 1.0796785   0.6620937 1.63070348 0.10295290      EBV 1.0796785 -1.005123
    ## 5 0.7279782   0.7319215 0.99461242 0.31992481  sex_fem 0.7279782 -4.594646
    ## 6 0.3816576   0.5335581 0.71530641 0.47441973     PDL1 0.3816576  2.616437
    ##            1 MeanDecreaseAccuracy MeanDecreaseGini
    ## 1  0.0000000            0.0000000        0.5884010
    ## 2  8.4643898            8.9514974        1.6613782
    ## 3  0.4536143            1.7002216        3.5860663
    ## 4  1.0805885            0.4480967        0.8550596
    ## 5  0.9889172           -2.9122258        0.5571359
    ## 6 -0.3073140            1.7497756        2.4945216

``` r
# ---- SELECT SORTING KEY ----
# Default
sorting_key <- merged_df$abs_coef
# sorting_key <- merged_df$MeanDecreaseAccuracy
# sorting_key <- merged_df$MeanDecreaseGini

# ---- Sort ----
merged_df$variable <- factor(merged_df$variable, levels = merged_df$variable[order(sorting_key, decreasing = FALSE)])

# ---- Adjust legend's location in Logit.Reg. ----
legend_pos <- c(1.2, 0.2)  # ex) "bottom", "left", c(0.9, 0.1)

# ---- Plot: Logistic Regression ----
p_logit <- ggplot(merged_df, aes(x = variable, y = Estimate)) +
  geom_col(aes(fill = Estimate > 0), width = 0.6) +
  coord_flip() +
  # geom_text(aes(label = sprintf("%.2f", Estimate)), 
  #           vjust = ifelse(merged_df$Estimate > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Positively\ncorrelated to\nReversePD == 1")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = legend_pos,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Coefficient",
       y = "Estimate")

# ---- Plot: RF Accuracy (Modified to resemble varImpPlot) ----
p_mda <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseAccuracy), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Accuracy",
       y = "MeanDecreaseAccuracy")

# ---- Plot: RF Gini (Modified to resemble varImpPlot) ----
p_gini <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseGini)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseGini), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Gini",
       y = "MeanDecreaseGini")


# ---- Combine all ----
library(grid)
library(gridExtra)

# Text labels
title_all <- textGrob("Scaled Variables' indices related to Target", gp = gpar(fontsize = 14, fontface = "bold"),
                      just = 'center') # center-align in entire plot
predictor_var_text <- textGrob("Target: ReversePD(=1; PD=0)", gp = gpar(fontsize = 12),
                               just = 'center')
title_logit <- textGrob("Logistic Regression", gp = gpar(fontsize = 13, fontface = "bold"))
title_rf    <- textGrob("Random Forest", gp = gpar(fontsize = 13, fontface = "bold"))



# subtitle row: 1st column = title_logit, 2nd+3rd column merged = title_rf
subtitle_row <- arrangeGrob(
  grobs = list(title_logit, title_rf),
  ncol = 3,
  layout_matrix = rbind(c(1, 2, 2))
)

# Plot row: actual plots
plot_row <- arrangeGrob(
  grobs = list(p_logit, p_mda, p_gini),
  ncol = 3,
  widths = c(2.5, 2.5, 2.5)
)

# Combine titles and plots
grid.arrange(
  title_all,
  predictor_var_text,
  subtitle_row,
  plot_row,
  nrow = 4,
  heights = c(0.12, 0.08, 0.12, 1)  # control height of title vs plot area
)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-127-1.png)<!-- -->

### Target: DCB

#### Logit.Reg.+LASSO

``` r
# Apply logistic regression with LASSO for 'Death' prediction
library(glmnet)

# Prepare predictor matrix and response vector
x <- as.matrix(df_surv_alive_included[, !names(df_surv_alive_included) %in% c("Responder", "PD", "PFS", "DiseaseDiagnosed", "Death", "Alive", "OS", "HER2",
                                                                "DCB")])    # Might be weird to use post-knowledge
y <- as.numeric(df_surv_alive_included$DCB)

# alpha = 1: LASSO, 0.5: Elastic Net
fit <- glmnet(x, y, family = "gaussian", # gaussian means logistic regression on regression problem
              type.measure = "mse",
              standardize = TRUE, # standardize the data -> Can differ the coefficient's axis
              # standardize = FALSE,
              alpha = 1)  # alpha=1 means LASSO, 0.5 means Elastic Net

# 계수 추출 (glmnet predict를 이용)
beta_pred <- predict(fit, type = "coefficients")
lambda_seq <- fit$lambda

# 계수 tidy하게 변환
coef_list <- lapply(1:length(lambda_seq), function(i) {
  beta <- as.matrix(beta_pred[, i])  # i번째 lambda에 대한 계수
  df <- data.frame(variable = rownames(beta),
                   coefficient = as.numeric(beta),
                   lambda = lambda_seq[i])
  df
})

coef_long <- bind_rows(coef_list) |>
  filter(variable != "(Intercept)")  # Exclude the intercept

ggplot(coef_long, aes(x = log(lambda), y = coefficient, color = variable)) +
  geom_line(linewidth = 1.2) +  # Increased line thickness
  theme_minimal() +
  labs(
      # title = "Elastic Net Regularization Path",
      title = "LASSO Regularization Path of predicing DCB",
       x = "log(lambda)",
       y = "Coefficient") +
  theme(legend.position = "right")
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-128-1.png)<!-- -->
\#### Only logit (No LASSO)

``` r
# Fit logistic reg. with scaled x
glm_fit <- glm(y ~ scale(x), family = binomial())

coef_summary <- summary(glm_fit)$coefficients
rownames(coef_summary) <- gsub("scale(x)", "", rownames(coef_summary), fixed = TRUE)
coef_summary
```

    ##               Estimate  Std. Error      z value  Pr(>|z|)
    ## (Intercept) -0.1972551  92.2917525 -0.002137299 0.9982947
    ## EBV          1.5034318   0.7458331  2.015775129 0.0438235
    ## MMRD         0.6840290   0.6018301  1.136581569 0.2557132
    ## age          1.1998524   0.8493548  1.412663339 0.1577547
    ## ATM_mut      3.5296133 470.5895000  0.007500408 0.9940156
    ## PDL1         0.2484918   0.5650450  0.439773503 0.6601012
    ## sex_fem      0.7735084   0.7926325  0.975872619 0.3291276

``` r
coef_df <- as.data.frame(coef_summary)
coef_df$variable <- rownames(coef_df)
# Store rows except of 'intercept'
coef_df <- coef_df[coef_df$variable != "(Intercept)", ]

# Sort variables by descending order of their absolutes
coef_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  arrange(desc(abs_coef))

ggplot(coef_df, aes(x = reorder(variable, abs_coef), y = Estimate)) +
  geom_col(aes(fill = Estimate > 0)) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance (Logistic Regression Coefficients)",
       x = "Variable",
       y = "Coefficient") +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Direction"))
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-130-1.png)<!-- -->
\#### Random Forest

``` r
library(randomForest)
# Fit a Random Forest model
rf <- randomForest(y = factor(y), x = scale(x), importance = TRUE)

# Display variable importance
importance(rf)
```

    ##                  0          1 MeanDecreaseAccuracy MeanDecreaseGini
    ## EBV      5.7217881  7.3105839             7.803361        1.7248729
    ## MMRD    -6.1833988 -2.9718720            -6.189024        0.4783830
    ## age      2.7558814  3.9779348             3.946534        3.7020150
    ## ATM_mut  0.0000000  0.0000000             0.000000        0.6010081
    ## PDL1    -0.5895073 -3.5274437            -2.149183        2.3324222
    ## sex_fem -4.0760319 -0.7245857            -3.706063        0.5141243

``` r
varImpPlot(rf)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-132-1.png)<!-- -->
\#### Combine LASSO with Logit&RF

``` r
# Extract importance measures
rf_importance <- as.data.frame(importance(rf))
rf_importance$variable <- rownames(rf_importance)

# Merge all information by variable
merged_df <- coef_df |>
  mutate(abs_coef = abs(Estimate)) |>
  left_join(rf_importance, by = "variable")
merged_df
```

    ##    Estimate  Std. Error     z value  Pr(>|z|) variable  abs_coef          0
    ## 1 3.5296133 470.5895000 0.007500408 0.9940156  ATM_mut 3.5296133  0.0000000
    ## 2 1.5034318   0.7458331 2.015775129 0.0438235      EBV 1.5034318  5.7217881
    ## 3 1.1998524   0.8493548 1.412663339 0.1577547      age 1.1998524  2.7558814
    ## 4 0.7735084   0.7926325 0.975872619 0.3291276  sex_fem 0.7735084 -4.0760319
    ## 5 0.6840290   0.6018301 1.136581569 0.2557132     MMRD 0.6840290 -6.1833988
    ## 6 0.2484918   0.5650450 0.439773503 0.6601012     PDL1 0.2484918 -0.5895073
    ##            1 MeanDecreaseAccuracy MeanDecreaseGini
    ## 1  0.0000000             0.000000        0.6010081
    ## 2  7.3105839             7.803361        1.7248729
    ## 3  3.9779348             3.946534        3.7020150
    ## 4 -0.7245857            -3.706063        0.5141243
    ## 5 -2.9718720            -6.189024        0.4783830
    ## 6 -3.5274437            -2.149183        2.3324222

``` r
# ---- SELECT SORTING KEY ----
# Default
sorting_key <- merged_df$abs_coef
# sorting_key <- merged_df$MeanDecreaseAccuracy
# sorting_key <- merged_df$MeanDecreaseGini

# ---- Sort ----
merged_df$variable <- factor(merged_df$variable, levels = merged_df$variable[order(sorting_key, decreasing = FALSE)])

# ---- Adjust legend's location in Logit.Reg. ----
legend_pos <- c(1, 0.3)  # ex) "bottom", "left", c(0.9, 0.1)

# ---- Plot: Logistic Regression ----
p_logit <- ggplot(merged_df, aes(x = variable, y = Estimate)) +
  geom_col(aes(fill = Estimate > 0), width = 0.6) +
  coord_flip() +
  # geom_text(aes(label = sprintf("%.2f", Estimate)), 
  #           vjust = ifelse(merged_df$Estimate > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"),
                    guide = guide_legend(title = "Positively\ncorrelated to\nDCB == 1")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = legend_pos,
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Coefficient",
       y = "Estimate")

# ---- Plot: RF Accuracy (Modified to resemble varImpPlot) ----
p_mda <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseAccuracy)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseAccuracy), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Accuracy",
       y = "MeanDecreaseAccuracy")

# ---- Plot: RF Gini (Modified to resemble varImpPlot) ----
p_gini <- ggplot(merged_df, aes(x = variable, y = MeanDecreaseGini)) +
  geom_segment(aes(xend = variable, y = 0, yend = MeanDecreaseGini), 
               linetype = "dashed", color = "gray30", size = 1) +
  geom_point(color = "black", size = 3.5) +
  coord_flip() +
  theme_minimal(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mean Decrease Gini",
       y = "MeanDecreaseGini")


# ---- Combine all ----
library(grid)
library(gridExtra)

# Text labels
title_all <- textGrob("Scaled Variables' indices related to Target", gp = gpar(fontsize = 14, fontface = "bold"),
                      just = 'center') # center-align in entire plot
predictor_var_text <- textGrob("Target: DCB(=1; NCB=0)", gp = gpar(fontsize = 12),
                               just = 'center')
title_logit <- textGrob("Logistic Regression", gp = gpar(fontsize = 13, fontface = "bold"))
title_rf    <- textGrob("Random Forest", gp = gpar(fontsize = 13, fontface = "bold"))



# subtitle row: 1st column = title_logit, 2nd+3rd column merged = title_rf
subtitle_row <- arrangeGrob(
  grobs = list(title_logit, title_rf),
  ncol = 3,
  layout_matrix = rbind(c(1, 2, 2))
)

# Plot row: actual plots
plot_row <- arrangeGrob(
  grobs = list(p_logit, p_mda, p_gini),
  ncol = 3,
  widths = c(2.5, 2.5, 2.5)
)

# Combine titles and plots
grid.arrange(
  title_all,
  predictor_var_text,
  subtitle_row,
  plot_row,
  nrow = 4,
  heights = c(0.12, 0.08, 0.12, 1)  # control height of title vs plot area
)
```

![](%5BVT%5Dfind_genes_related_to_long_surv_files/figure-gfm/unnamed-chunk-134-1.png)<!-- -->
