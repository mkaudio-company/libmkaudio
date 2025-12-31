#include <cmath>
#include <algorithm>

#include "saturation.h"

namespace MKAudio
{
Saturation::Saturation(double alphaPlus, double alphaMinus, 
                   double betaPlus, double betaMinus, 
                   double deltaBias, bool flip) 
{
    // 1. Set Alpha (Knee/Drive)
    driveAlphaPlus = std::max(alphaPlus, 1e-4);
    driveAlphaMinus = std::max(alphaMinus, 1e-4);
 
    // 2. Set Beta (Compression/Gain)
    compressionBetaPlus = betaPlus;
    compressionBetaMinus = betaMinus;

    // 3. Set Bias (Delta)
    biasDelta = deltaBias;

    // 4. Set Polarity (Gamma - Boolean)
    flipPolarity = flip;

    // 5. Recalculate Normalization Factors
    normFactorPlus = 1.0 / std::log2(1.0 + driveAlphaPlus);
    normFactorMinus = 1.0 / std::log2(1.0 + driveAlphaMinus);
}

double Saturation::process(double inputSample)
{
    
    double outputValue;

    // 1. Calculate Saturated Output (S(x, δ))
    if (inputSample >= biasDelta)
    {
        // Processing Positive side (x >= δ)
        double relativeInput = inputSample - biasDelta; 
        double logOut = std::log2(1.0 + (driveAlphaPlus * relativeInput));
        outputValue = compressionBetaPlus * (logOut * normFactorPlus); 
    }
    else
    {
        // Processing Negative side (x < δ)
        double relativeInput = biasDelta - inputSample; 
        double logOut = std::log2(1.0 + (driveAlphaMinus * relativeInput));
        // Output must be negative before final flip
        outputValue = -1.0 * compressionBetaMinus * (logOut * normFactorMinus);
    }
    
    // 2. Apply Boolean Polarity Flip (γ)
    if (flipPolarity) return -outputValue;
    return outputValue;
}

} // namespace MKAudio