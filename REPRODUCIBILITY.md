# Reproducibility

This repository is organized into loader scripts, run scripts, and visualization scripts to enable reproducibility of all analyses presented in the HEAP manuscript.

## Loader Scripts

The following scripts load and structure data into custom S4 objects required for analysis:

- `loaderProt_PGS_PXSv2.R`  
- `loaderProt_PGS_PXS_mediationv2.R`

These scripts construct two primary objects:

### `PXSconstruct`

```r
PXSconstruct <- setClass(
  "PXSconstruct",
  slots = c(
    Elist = "list",           # Exposure dataframes separated by category
    Elist_names = "character",# Names of each exposure category
    Eid_cat = "data.frame",   # Mapping of exposure IDs to categories
    ordinalIDs = "character", # List of ordinal exposure variable names

    UKBprot_df = "data.frame",# Protein variable dataframe
    protIDs = "character",    # Protein variable names

    covars_df = "data.frame", # Covariate dataframe
    covars_list = "character" # Covariate variable names
  )
)
```

### `PXSmediation`

```r
PXSmediation <- setClass(
  "PXSmediation",
  slots = c(
    Protlist = "list",        # Protein dataframes by category
    protIDs = "character",    # Protein names

    DZ_df = "data.frame",     # Disease-related data (first occurrence, age, death, etc.)
    DZ_ids = "character",     # Disease identifiers for age of onset

    covars_df = "data.frame", # Covariate dataframe
    covars_list = "character" # Covariate variable names
  )
)
```

To apply HEAP to any dataset, the relevant data must be loaded into the appropriate slots of these structures.

## Run Scripts

The following scripts perform core HEAP analyses:

- `runProt_PXS_PXS_finalv3.R`: Runs the variance decomposition module.  
- `runProt_PGS_univariae_v3.R` and `runProt_PGS_univarInteract_v2.R`: Run the gene-by-environment association module.  
- `runProt_PGS_PXS_mediationv4.R`: Runs the mediation analysis module.  
- `mediation_bootstrap_run.R` and `mediation_bootstrap_source.R`: Support bootstrapping for mediation analysis.

## Visualization Scripts

Visualization scripts are located in the `Visualizations/` folder and follow a modular design:

**Step 1: Build Data Structures**

- `analyHEAP_association.R`  
- `analyHEAP_intervention.R`  
- `analyHEAP_mediation.R`  
- `analyProt_PGS_PXS_final.R`  

These scripts generate the intermediate objects required for visualization.

**Step 2: Create Figures**

Visualizations used to generate figures in the HEAP manuscript are located in the `Module1` to `Module4` and `ModuleExt` folders. Each folder is organized by topic and type of plot (e.g., scatter plots, bar charts, interactive plots, tables), and contains its own `README` describing the scripts inside.

## Contact

For further information or questions, please contact: shakson_isaac@g.harvard.edu
