THIS FOLDER CONTAINS THE CODE, written in Matlab (version R2021a) and Jupyter Notebooks.

To run the scripts, download the folder structure as given here with */code/, */var/, and */data/ sub-folders, 
and set the path (and parameters, if you want) at the beginning of the scripts. 

IF YOU HAVE ANY QUESTIONS, DO NOT HESITATE TO ASK:

Isa Dallmer-Zerbe
Institute of Computer Science
The Czech Academy of Sciences
Prague

dallmer-zerbe@cs.cas.cz


## Overview

This code implements a model-based approach to characterize different seizure onset patterns that:
1. Simulates EEG signals using the Wendling (2005) neural mass model across parameter sweeps
2. Extracts signal features from both simulated and real intracranial EEG data (seizure recordings from EPILEPSOAE database)
3. Classifies epileptic brain states (interictal, pre-onset, onset, ictal)
4. Estimates optimal model parameters by matching real EEG features to simulated ones
5. Compares identified model parameters across time (interictal to ictal) for different seizure onset patterns.

## Repository Structure

```
├── c1_get_features_parsweep.m    # Generate features from model parameter sweeps
├── c1_get_features_real.m        # Extract features from real EEG recordings
├── c2_get_model_result.m         # Classify brain states and fit ABG parameters
├── c2_prepare_outcome_vars.m     # Aggregate results across channels and time windows
├── c3_analyse_results.m          # Analysis and visualization scripts
├── c3_sz_types_classification.ipynb # Python notebook for seizure onset pattern classification
└── functions/                     # Helper functions 
    ├── Wendl2005stimnoise.m
    ├── calc_feat.m
    ├── generate_prototypes.m
    ├── classify_type.m
    └── optimize_abg.m
```

## Code Workflow

### Step 1: Generate Simulated Features (`c1_get_features_parsweep.m`)

Simulates the Wendling model across a parameter space and extracts 11 signal features:

**Parameters:**
- **A**: [2, 14], step 0.25
- **B**: [1, 25], step 0.5  
- **G**: [1, 30], step 0.5
- Simulation duration: 5 seconds
- Sampling rate: 512 Hz

**Features extracted (11 total):**
- Band power (5 bands): b0-b4
- Alpha difference: alphdiff
- Spike absolute: spikeabs
- Signal mean: sigmean
- Signal variance: sigvar
- Autocorrelation: autocorrel
- Line length: linelen

**Output:** `parsweep{N}.mat` containing z-normalized features for all (A,B,G) combinations

### Step 2: Extract Real EEG Features (`c1_get_features_real.m`)

Processes real seizure recordings from the EPILEPSIAE database:

**Processing:**
- Segments continuous EEG into 5-second windows
- Extracts the same 11 features per segment
- Z-normalizes features relative to interictal baseline (first 30 seconds of the recording)

**Input format:**
```matlab
% Seizure files: Seiz*.mat
data(channel_ind, time_ind)   % Multi-channel EEG
info                          % Metadata structure
```

**Output:** `FEATURES-epilepsiae.mat` structure array with fields:
- `fileID`: Seizure and channel identifier
- `featuresF0`: Normalized features (segments × 11)
- `featuresF0_raw`: Raw feature values

### Step 3: Brain State Classification & Parameter Fitting (`c2_get_model_result.m`)

Performs the core analysis by comparing real and simulated features:

**Key processes:**
1. **Prototype generation**: Creates cluster centroids for brain state classification
2. **PCA projection**: Reduces feature dimensionality (default: 4 components)
3. **Brain state classification**: Assigns each segment to one of 4 brain states
4. **ABG optimization**: Finds parameter set with minimum feature distance and thus "optimal ABG" (indicative of excitation-inhibition levels)

**Configuration:**
```matlab
pars.nclust = 4;           % Number of brain state clusters
pars.pcafeat = true;       % Use PCA preprocessing
pars.n_comp = 4;           % PCA components
pars.parsweepnum = 2;      % Which parameter sweep to use
```

**Output:** `RESULT-epilepsiae.mat` structure with:
- `minerrtype`: Classified brain state per segment
- `minerrparsABG`: Optimal (A, B, G) parameters per segment (indicative of excitation-inhibition levels)
- `minerrfeatures`: Simulated features at optimal parameters

### Step 4: Aggregate Outcome Variables (`c2_prepare_outcome_vars.m`)

Summarizes results across seizure phases and brain regions:

**Time windows:**
- Interictal: 60-30s before seizure onset
- Pre-onset: 30-0s before onset
- Onset: 0-10s after onset
- Ictal: 10-25s after onset

**Spatial aggregation:**
- All channels
- SOZ (seizure onset zone) channels only
- Non-SOZ channels

**Output:** `OutcomeVars-epilepsiae.mat` with averaged A, B, G, and brain state values

## Requirements

### MATLAB Dependencies
- MATLAB R2019b or later
- Parallel Computing Toolbox (for `parfor` loops)
- Signal Processing Toolbox

### Python Dependencies (for notebook)
```bash
pip install numpy pandas matplotlib scikit-learn jupyter
```

### External Functions
The pipeline requires custom helper functions in the `functions/` directory:
- `Wendl2005stimnoise.m`: Wendling model simulation
- `calc_feat.m`: Feature extraction
- `generate_prototypes.m`: Cluster prototype generation
- `classify_type.m`: Brain state classification
- `optimize_abg.m`: Parameter optimization
- `get_types.m`: Time window extraction
- `type_labeling.m`: Helper function for brain state classification
- `preprocess_features.m`: Helper function for cluster prototype generation

## Installation & Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/epilepsy-analysis.git
cd epilepsy-analysis
```

2. Set up paths in each script:
```matlab
PATH = 'C:\your\project\path';
PATH_vars = [PATH, '\vars'];
PATH_code = [PATH, '\functions'];
PATH_data = [PATH, '\data'];
```

3. Add helper functions to `functions/` directory

4. Prepare data in expected format (see Data Format section: Input)

## Data Format

### Input: Seizure Data
```matlab
% File: Seiz{seizID}_{patientID}.mat
data        % [channels × timepoints] EEG array
info.fs     % Sampling frequency
info.SOZ    % [1 × channels] logical, seizure onset zone markers
info.seiz_start_index  % Seizure onset sample
```

### Output: Feature Structure
```matlab
FEATURES(file_ind, channel_ind)
    .fileID           % Identifier string
    .featuresF0       % [segments × 11] normalized features
    .featuresF0_raw   % [segments × 11] raw features
```

### Output: Result Structure
```matlab
RESULT(file_ind, channel_ind)
    .minerrtype       % [segments × 1] brain state (1-4)
    .minerrparsABG    % [segments × 3] fitted (A, B, G)
    .minerr           % [segments × 1] fitting error
```

## Usage Example

```matlab
%% 1. Generate parameter sweep features
cd('/path/to/project')
c1_get_features_parsweep  % Creates parsweep2.mat

%% 2. Extract features from real data
c1_get_features_real      % Creates FEATURES-epilepsiae.mat

%% 3. Classify and fit parameters
c2_get_model_result       % Creates RESULT-epilepsiae.mat

%% 4. Aggregate outcomes
c2_prepare_outcome_vars   % Creates OutcomeVars-epilepsiae.mat

%% 5. Analyze results
c3_analyze_results        % Plots main figures of the manuscripts and 

```

## Parallel Processing Configuration

The code uses MATLAB's parallel processing for computational efficiency:

```matlab
% In c1_get_features_parsweep.m
parpool;  % Default workers

% In c2_get_model_result.m
parpool('local', 2);  % Customize worker count
```

Adjust the number of workers based on your system resources.

## Performance Considerations

- **Parameter sweep**: ~50,000 simulations for default grid (parallelized)
- **Feature extraction**: ~1-2 seconds per channel per seizure
- **Model fitting**: ~5-10 seconds per channel (depends on sweep size)
- **Memory**: Convert to `single` precision for large datasets

## Citation

If you use this code, please cite:

```
Dallmer-Zerbe, I., et al. (2023). Model-based brain state estimation for epilepsy.
[Publication details]
```

## References

- Wendling, F., et al. (2005). Epileptic fast activity can be explained by a model of impaired GABAergic dendritic inhibition. *European Journal of Neuroscience*.
- Fietkiewicz, C., & Loparo, K. A. (2016). Stochastic neural field model of stimulus-dependent variability in cortical neurons.

## Acknowledgments

This work uses data from the EPILEPSIAE database. Due to data restrictions, no data or output structures are provided. 
