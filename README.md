# **Stimulus-Response Model in MATLAB**

## **Overview**  
This repository contains MATLAB code implementing a stimulus-response model based on the linear–nonlinear-linear–nonlinear (LN-LN) cascade model (referencing: Heeger et al., 2015). The model is designed to simulate the neural population's response to various Gabor stimuli, commonly used in visual neuroscience.

The primary goal of this model is to investigate trial-to-trial variability in neural responses to identical stimuli. Specifically, it explores this variability through the framework of choice probability (CP), which measures the correlation between neural activity and behavioral decisions. By simulating and analyzing these responses, the model allows us to:

- Understand how CP is influenced by various factors, such as stimulus properties, behavioral biases, etc.
- Explore how these factors affect the relationship between sensory inputs and decision-making processes.
- Identify potential confounding influences that may skew CP measurements, which is crucial for improving the design of behavioral and neural experiments.

Through this analysis, the project aims to provide insights into refining experimental setups to more accurately capture neural activity related to decision-making, and to better handle sources of noise or bias that might impact CP.

### **Key Features**:
- Simulates neural responses to Gabor stimuli with varying orientations (adapted from LN-LN cascade model)
- Uses modulated poisson process to generate spikes across for neurons
- Implement a choice probability (CP) metric for analyzing trial-by-trial correlation between decision and neural firing
- Customizable gain profiles for individual neurons to mimic the top-down modulation of neural firing.

## **Installation**

1. **Clone the repository**:

    ```bash
    git clone https://github.com/avinash1598/StimRespModel.git
    ```

2. **Open the MATLAB code**:  
   Open MATLAB and navigate to the cloned repository folder.

3. **Run the script**:  
   The main script `GenerateTrialDataOptimized.m` implements the LN-LN (Linear-Nonlinear) cascade model and a modulated Poisson process to simulate neural responses to various stimulus orientations. The generated neural spike data is then saved to a file for further analysis.

## **Usage**

### **Main Scripts for Model**:

1. **`GenerateTrialDataOptimized.m`**:  
   Main script that implements the LN-LN cascade model and uses a modulated Poisson process to generate trial-by-trial neural response data. In this script, key simulation parameters, such as the number of neurons, stimulus vector, and stimulus duration, can be defined. The script then runs the simulation and generates the following key data, which is saved to a file:
   - trialMatrix: Behavioral data containing details of the stimulus orientations used for the simulation and the decision made by the model for each of those orientations.
   - trialResponses: Neural responses, i.e., spiking data from all the defined neurons for each stimulus shown to the model.
   
2. **`GenerateTrialData.m`**:  
   Functions similarly to `GenerateTrialDataOptimized.m`, except that the code is not optimized for memory and might exhaust RAM as the simulation parameters increase. Whereas, `GenerateTrialDataOptimized.m` is optimized for memory and can run on systems with very low RAM (8 GB).
   
3. **`decodeOrientationFromSpikes.m`**:  
   - TBD

4. **`generateModulatedPoissonSpikes.m`**:  
   - TBD
    
5. **`orientationTunedFiringRate.m`**:  
   - TBD
  
6. **`getGainProfile.m`**:
   - TBD
