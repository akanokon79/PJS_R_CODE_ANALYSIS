 PJS_BFLQR_ANALYSIS

R code and data for the simulation study and real-data analysis accompanying a manuscript submitted to the Pakistan Journal of Statistics.

# Files

PJS_Simulation_Code.R — Simulation study.
PJS_Real_Data_Analysis.R — Real-data analysis.
thesisdata.csv — Dataset used for the real-data analysis.

# Requirements

The real-data script uses the R packages `quantreg`, `MASS`, `GIGrvg`, and `coda`.

Install them using:
install.packages(c("quantreg", "MASS", "GIGrvg", "coda"))

# How to Run

Place all three files in the same folder, open R or RStudio in that folder, and run:
source("PJS_Simulation_Code.R")
source("PJS_Real_Data_Analysis.R")

The scripts use the supplied dataset and define their respective analysis settings.

# Reproducibility

The code and dataset are provided to support reproducibility of the reported analyses. Results should be verified by running the scripts with the required packages and computational settings.

# Data Description

The file thesisdata.csv contains the response variable and ten predictors used in the real-data analysis.

Y — CO₂ emissions (response variable).
X1 — GDP per capita.
X2 — GDP growth.
X3 — Industry value added.
X4 — Foreign direct investment (FDI) inflows.
X5 — Natural resources rents.
X6 — Energy use.
X7 — Electric power consumption.
X8 — Access to electricity.
X9 — Renewable energy consumption.
X10 — Urban population growth
Retain the variable names Y and X1–X10 in the CSV file to ensure compatibility with the analysis script.
