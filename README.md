# **Stimulus-Response Model in MATLAB**

## **Overview**  
This repository contains MATLAB code implementing a stimulus-response model based on the linear–nonlinear-linear–nonlinear (LN-LN) cascade model (referencing: Heeger et al., 2015). The model is designed to simulate the neural population's response to various Gabor stimuli, commonly used in visual neuroscience.

The primary goal of this model is to investigate trial-to-trial variability in neural responses to identical stimuli. Specifically, it explores this variability through the framework of choice probability (CP), which measures the correlation between neural activity and behavioral decisions. By simulating and analyzing these responses, the model allows us to:

- Understand how CP is influenced by various factors, such as stimulus properties, behavioral biases, etc.
- Explore how these factors affect the relationship between sensory inputs and decision-making processes.
- Identify potential confounding influences that may skew CP measurements, which is crucial for improving the design of behavioral and neural experiments.

Through this analysis, the project aims to provide insights into refining experimental setups to more accurately capture neural activity related to decision-making, and to better handle sources of noise or bias that might impact CP.

### **Key Features**:
- Simulates neural responses to Gabor stimuli with varying orientations.
- Computes firing rates across neurons over time.
- Implements a choice probability (CP) metric for decision-related analysis.
- Supports statistical significance testing via permutation tests.
- Customizable gain profiles for individual neurons.

## **Installation**

1. **Clone the repository**:

    ```bash
    git clone https://github.com/your-username/stimulus-response-model.git
    ```

2. **Open the MATLAB code**:  
   Open MATLAB and navigate to the cloned repository folder.

3. **Run the script**:  
   You can run the main script `main.m` or any other provided functions directly in MATLAB.

## **Usage**

### **Main Scripts and Functions**:

1. **`main.m`**:  
   - The entry point to the model simulation. It sets up the stimulus parameters and executes the analysis.
   
2. **`getGainProfile.m`**:  
   - Generates gain profiles for a given number of neurons over time using cumulative Gaussian distributions.
   
3. **`computeChoiceProbability.m`**:  
   - Calculates the choice probability (CP) for each neuron based on its firing rates and trial labels.

4. **`permutationTest.m`**:  
   - Performs a permutation test to assess the significance of the computed CP values.

### **Example Usage**:

```matlab
nNeurons = 100;
timeBins = linspace(0, 1, 100); % 100 time bins for 1-second stimulus duration
stimDur = 1; % Stimulus duration in seconds

% Generate gain profiles
gainProfiles = getGainProfile(nNeurons, timeBins, stimDur);

% Compute choice probability for simulated data
cpValues = computeChoiceProbability(firingRates, labels);

% Assess significance using permutation test
pValues = permutationTest(firingRates, labels, numPermutations);
