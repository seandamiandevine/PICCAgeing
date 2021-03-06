# PICCAgeing
Devine, Neumann, Levari, Wilson, &amp; Eppinger (2021), prevalence-induced concept change in human aging. Task code and data backup. 

## Directories
``` 
.
.
├── Analysis.R
├── data
│   ├── Dots
│   └── Ethics
├── Models
│   ├── DDM
│   └── seq
├── PICCOA_2021.RData
├── README.md
└── TaskCode
    ├── EthicPICC.py
    ├── PICCcolour.py
    └── stimuli


```

### **data**
Contains subject data. Id # < 100 is OA, > 100 is YA. Dots and Ethics task data are saved in their own directories (as .csv)

### **Models**
Contains model code and individual model fits for both the DDM (/DDM/, Ratcliff & McKoon, 2008) and sequential decision-making model (/seq/, Wilson, 2018). 

### **TaskCode**
PsychoPy scrips for the the experiment. PICCcolour.py is the Dots Task. EthicPICC.py is the Ethics task. /stimuli/ contains all the stimuli used in the Ethics task--i.e., the research scenarios. 

## Standalone Files
* Analysis.R is the main analysis script
* PICCOA_2021.RData is the R environment for an executed Analysis.R
* README.md is this file


