#pragma once

namespace MKAudio
{

class Saturation
{
private:
    // 1. Asymmetry/Knee parameters (Alpha)
    double driveAlphaPlus, driveAlphaMinus;  
    // 2. Compression/Gain parameters (Beta)
    double compressionBetaPlus, compressionBetaMinus;
    // 3. Bias parameter (Delta)
    double biasDelta; 
    
    // 4. Polarity flip parameter (Gamma, now Boolean)
    bool flipPolarity; // true = invert, false = no inversion
    
    // Normalization factors (Calculated once)
    double normFactorPlus, normFactorMinus;

public:
    Saturation(double alphaPlus, double alphaMinus, 
               double betaPlus, double betaMinus, 
               double deltaBias, bool flip);

    // The Final Core Asymmetrical, Biased "Numeric Modeling" Algorithm
    double process(double inputSample);
};

} // namespace MKAudio