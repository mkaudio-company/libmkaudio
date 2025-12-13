#pragma once

namespace MKAudio
{

class BiquadFilter
{
protected:
    // Filter Coefficients
    double b0 = 1.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0; // Note: a0 is always 1.0, so we normalize by it

    // State Variables (delay elements)
    double x1 = 0.0, x2 = 0.0; // Input delays: x[n-1], x[n-2]
    double y1 = 0.0, y2 = 0.0; // Output delays: y[n-1], y[n-2]

    double sampleRate;

public:
    // Constructor: Takes the audio sampling rate (e.g., 44100.0)
    BiquadFilter(double fs) : sampleRate(fs) {}

    // Pure virtual function: Forces derived classes to implement
    // the logic for calculating coefficients.
    virtual void calculateCoefficients(double fc, double Q, double gain_dB) = 0;

    // The core DSP processing function
    double process(double input);
};

// --- 1. Low-Pass Filter (LP) ---
// Second-order (12 dB/octave) Low-Pass
class LowPassFilter : public BiquadFilter {
public:
    LowPassFilter(double fs) : BiquadFilter(fs) {}

    // Q is the resonance factor (e.g., 0.707 for Butterworth)
    void calculateCoefficients(double fc, double Q, double /* gain_dB unused */) override;
};

// --- 2. High-Pass Filter (HP) ---
// Second-order (12 dB/octave) High-Pass
class HighPassFilter : public BiquadFilter {
public:
    HighPassFilter(double fs) : BiquadFilter(fs) {}

    void calculateCoefficients(double fc, double Q, double /* gain_dB unused */) override;
};

// --- 3. Bell (Peaking) Filter ---
class BellFilter : public BiquadFilter {
public:
    BellFilter(double fs) : BiquadFilter(fs) {}

    // Q controls the bandwidth; gain_dB controls the boost/cut
    void calculateCoefficients(double fc, double Q, double gain_dB) override;
};

// --- 4. Low Shelf Filter ---
class LowShelfFilter : public BiquadFilter {
public:
    LowShelfFilter(double fs) : BiquadFilter(fs) {}

    // Q is often fixed to 0.707 but can be used to control the transition slope
    void calculateCoefficients(double fc, double Q, double gain_dB) override;
};

// --- 5. High Shelf Filter ---
class HighShelfFilter : public BiquadFilter {
public:
    HighShelfFilter(double fs) : BiquadFilter(fs) {}

    void calculateCoefficients(double fc, double Q, double gain_dB) override;
};

} // namespace MKAudio