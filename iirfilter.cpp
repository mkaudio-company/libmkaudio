#include <cmath>
#include <algorithm>

#include "iirfilter.h"

namespace MKAudio
{

double BiquadFilter::process(double input)
{
    // The Biquad difference equation:
    // y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
    double output = b0 * input + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;

    // Update state variables for the next sample
    x2 = x1;
    x1 = input;
    y2 = y1;
    y1 = output;

    return output;
}

void LowPassFilter::calculateCoefficients(double fc, double Q, double /* gain_dB unused */)
{
    // Clamp frequency and Q for stability
    fc = std::min(fc, sampleRate / 2.0 - 1.0);
    Q = std::max(0.5, Q); 

    // 1. Pre-warping and intermediate variables
    double w0 = 2.0 * M_PI * fc / sampleRate;
    double alpha = std::sin(w0) / (2.0 * Q);
    double cosw0 = std::cos(w0);

    // 2. Denominator normalization factor
    double a0_inv = 1.0 / (1.0 + alpha);

    // 3. Coefficients calculation (LPF transfer function)
    b0 = ( (1.0 - cosw0) / 2.0 ) * a0_inv;
    b1 = ( 1.0 - cosw0 ) * a0_inv;
    b2 = ( (1.0 - cosw0) / 2.0 ) * a0_inv;
    a1 = ( -2.0 * cosw0 ) * a0_inv;
    a2 = ( 1.0 - alpha ) * a0_inv;
}

void HighPassFilter::calculateCoefficients(double fc, double Q, double /* gain_dB unused */)
{
    // Clamp frequency and Q for stability
    fc = std::min(fc, sampleRate / 2.0 - 1.0);
    Q = std::max(0.5, Q); 
    
    // 1. Pre-warping and intermediate variables
    double w0 = 2.0 * M_PI * fc / sampleRate;
    double alpha = std::sin(w0) / (2.0 * Q);
    double cosw0 = std::cos(w0);

    // 2. Denominator normalization factor
    double a0_inv = 1.0 / (1.0 + alpha);

    // 3. Coefficients calculation (HPF transfer function)
    b0 = ( (1.0 + cosw0) / 2.0 ) * a0_inv;
    b1 = ( -(1.0 + cosw0) ) * a0_inv;
    b2 = ( (1.0 + cosw0) / 2.0 ) * a0_inv;
    a1 = ( -2.0 * cosw0 ) * a0_inv;
    a2 = ( 1.0 - alpha ) * a0_inv;
}

void BellFilter::calculateCoefficients(double fc, double Q, double gain_dB)
{
    // Clamp frequency and Q for stability
    fc = std::min(fc, sampleRate / 2.0 - 1.0);
    Q = std::max(0.5, Q); 

    // 1. Pre-warping and intermediate variables
    double w0 = 2.0 * M_PI * fc / sampleRate;
    double A = std::pow(10.0, gain_dB / 40.0); // Linear gain for the boost/cut (sqrt(10^(dB/20)))
    double alpha = std::sin(w0) / (2.0 * Q);
    double cosw0 = std::cos(w0);

    // 2. Denominator normalization factor
    double a0_inv = 1.0 / (1.0 + alpha);

    // 3. Coefficients calculation (Bell transfer function)
    b0 = ( 1.0 + alpha * A ) * a0_inv;
    b1 = ( -2.0 * cosw0 ) * a0_inv;
    b2 = ( 1.0 - alpha * A ) * a0_inv;
    a1 = ( -2.0 * cosw0 ) * a0_inv;
    a2 = ( 1.0 - alpha ) * a0_inv;
}

void LowShelfFilter::calculateCoefficients(double fc, double Q, double gain_dB)
{
    fc = std::min(fc, sampleRate / 2.0 - 1.0);
    Q = std::max(0.5, Q); // Q controls the transition slope/peaking at fc

    // 1. Pre-warping and intermediate variables
    double w0 = 2.0 * M_PI * fc / sampleRate;
    double A = std::pow(10.0, gain_dB / 40.0);
    double cosw0 = std::cos(w0);
    double sinw0 = std::sin(w0);
    double beta = std::sqrt(A) / Q; // Simplified formula for a smooth slope

    // 2. Denominator and Numerator components
    double Ap1 = A + 1.0;
    double Am1 = A - 1.0;

    // Simplified coefficient derivation for low-shelf
    double denom = Ap1 + Am1 * cosw0 + beta * sinw0;

    double b0_num = A * (Ap1 - Am1 * cosw0 + beta * sinw0);
    double b1_num = 2.0 * A * (Am1 - Ap1 * cosw0);
    double b2_num = A * (Ap1 - Am1 * cosw0 - beta * sinw0);

    double a1_num = -2.0 * (Am1 + Ap1 * cosw0);
    double a2_num = -(Ap1 + Am1 * cosw0 - beta * sinw0);

    // 3. Final Coefficients
    double a0_inv = 1.0 / denom;
    b0 = b0_num * a0_inv;
    b1 = b1_num * a0_inv;
    b2 = b2_num * a0_inv;
    a1 = a1_num * a0_inv;
    a2 = a2_num * a0_inv;
}

void HighShelfFilter::calculateCoefficients(double fc, double Q, double gain_dB)
{
    fc = std::min(fc, sampleRate / 2.0 - 1.0);
    Q = std::max(0.5, Q);

    // 1. Pre-warping and intermediate variables (same as Low Shelf)
    double w0 = 2.0 * M_PI * fc / sampleRate;
    double A = std::pow(10.0, gain_dB / 40.0);
    double cosw0 = std::cos(w0);
    double sinw0 = std::sin(w0);
    double beta = std::sqrt(A) / Q;

    // 2. Denominator and Numerator components
    double Ap1 = A + 1.0;
    double Am1 = A - 1.0;

    // Simplified coefficient derivation for high-shelf
    double denom = Ap1 - Am1 * cosw0 + beta * sinw0;

    double b0_num = A * (Ap1 + Am1 * cosw0 + beta * sinw0);
    double b1_num = -2.0 * A * (Am1 + Ap1 * cosw0);
    double b2_num = A * (Ap1 + Am1 * cosw0 - beta * sinw0);

    double a1_num = 2.0 * (Am1 - Ap1 * cosw0);
    double a2_num = -(Ap1 - Am1 * cosw0 - beta * sinw0);

    // 3. Final Coefficients
    double a0_inv = 1.0 / denom;
    b0 = b0_num * a0_inv;
    b1 = b1_num * a0_inv;
    b2 = b2_num * a0_inv;
    a1 = a1_num * a0_inv;
    a2 = a2_num * a0_inv;
}

} // namespace MKAudio