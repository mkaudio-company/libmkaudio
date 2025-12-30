#pragma once

#include <memory>
#include <array>
#include <algorithm>
#include <cmath>

namespace MKAudio
{

// ==========================================
// 1. Optimized Component Classes
// ==========================================
class Component {
protected:
    int nodeA, nodeB;

public:
    Component(int n1, int n2);
    virtual ~Component() = default;

    // 1. Preprocess Step: Contributes to the static Y Matrix
    virtual double getConductance(double dt) const = 0;

    // 2. Real-time Step: Contributes to the dynamic J Vector
    virtual double getCurrentSource(double dt) const = 0;

    // 3. Post-Process Step: Updates internal memory
    virtual void updateState(double vA, double vB, double dt) = 0;

    // 4. Check if component value can change dynamically
    virtual bool isDynamic() const;

    int getNodeA() const;
    int getNodeB() const;
};

class Resistor : public Component {
    double conductance;
public:
    Resistor(int n1, int n2, double r);

    double getConductance(double dt) const override;
    double getCurrentSource(double dt) const override;
    void updateState(double vA, double vB, double dt) override;

    void setResistance(double r);
    double getResistance() const;
};

class Capacitor : public Component {
    double capacitance;
    double prevVoltage;

public:
    Capacitor(int n1, int n2, double c);

    double getConductance(double dt) const override;
    double getCurrentSource(double dt) const override;
    void updateState(double vA, double vB, double dt) override;
};

class Inductor : public Component {
    double inductance;
    double prevCurrent;

public:
    Inductor(int n1, int n2, double l);

    double getConductance(double dt) const override;
    double getCurrentSource(double dt) const override;
    void updateState(double vA, double vB, double dt) override;
};

// ==========================================
// Potentiometer Class
// ==========================================
// A potentiometer is modeled as two resistors in series:
// nodeA --- R_upper --- wiper (nodeW) --- R_lower --- nodeB
// The wiper position (0.0 to 1.0) determines the ratio
// When position = 0.0: R_upper = totalR, R_lower = 0 (wiper at nodeB)
// When position = 1.0: R_upper = 0, R_lower = totalR (wiper at nodeA)
//
// For audio taper (logarithmic), use setPositionLog()
class Potentiometer : public Component {
private:
    int nodeWiper;
    double totalResistance;
    double position;        // 0.0 to 1.0
    double conductanceUpper;
    double conductanceLower;

    static constexpr double MIN_R = 1.0; // Minimum resistance to avoid division by zero

    void updateConductances();

public:
    // n1 = nodeA (top), n2 = nodeB (bottom), nW = wiper node
    Potentiometer(int n1, int n2, int nW, double totalR, double initialPosition = 0.5);

    // This component is dynamic (resistance changes with knob position)
    bool isDynamic() const override;

    // For a potentiometer, we need special handling - it's actually TWO resistors
    // We'll return the combined parallel conductance for stamping purposes
    // But the actual stamping needs to be done differently
    double getConductance(double dt) const override;

    double getCurrentSource(double dt) const override;
    void updateState(double vA, double vB, double dt) override;

    // Set wiper position (0.0 to 1.0, linear)
    void setPosition(double pos);

    // Set wiper position with audio/log taper
    // Attempt to model a typical audio taper pot
    void setPositionLog(double pos);

    // Set wiper position with reverse log taper
    void setPositionRevLog(double pos);

    double getPosition() const;
    int getWiperNode() const;
    double getConductanceUpper() const;
    double getConductanceLower() const;
    double getTotalResistance() const;
};

// ==========================================
// 2. Optimized Circuit Engine
// ==========================================
template <std::size_t NumDevices, std::size_t NumNodes>
class Circuit {
private:
    // Fixed Size Arrays for Optimization
    std::array<std::shared_ptr<Component>, NumDevices> components{};

    // Row-major flattened matrix for Y: Y[row][col] -> Y_static[row * NumNodes + col]
    std::array<double, NumNodes * NumNodes> Y_static{};
    std::array<double, NumNodes * NumNodes> Y_work{};   // Scratchpad for solver
    std::array<double, NumNodes> J{};                    // Current vector
    std::array<double, NumNodes> nodes{};                // Solution (Voltages)

    int currentDevices;
    double dt;

    // Internal Optimized Solver (Gaussian Elimination)
    // Solves Y * x = J, result stored in J
    void solveLinearSystem(int N);

public:
    Circuit(double sampleRate);

    void addComponent(std::shared_ptr<Component> c);

    // ==========================================
    // Preprocess: Builds the static Y Matrix
    // Call this ONCE before audio starts
    // ==========================================
    void preprocess(double impedence);

    // ==========================================
    // Real-time Audio Process
    // ==========================================
    double process(double inputVoltage, int probeNode);
};

// ==========================================
// 3. Dynamic Circuit Engine (supports real-time parameter changes)
// ==========================================
template <std::size_t NumDevices, std::size_t NumNodes>
class DynamicCircuit {
private:
    std::array<std::shared_ptr<Component>, NumDevices> components{};
    std::array<std::shared_ptr<Potentiometer>, NumDevices> potentiometers{}; // Track pots separately

    std::array<double, NumNodes * NumNodes> Y_work{};
    std::array<double, NumNodes> J{};
    std::array<double, NumNodes> nodes{};

    int currentDevices = 0;
    int currentPots = 0;
    double dt;
    double sourceImpedance = 10.0;

    void solveLinearSystem(int N);
    void stampComponent(std::shared_ptr<Component>& comp);
    void stampPotentiometer(std::shared_ptr<Potentiometer>& pot);

public:
    DynamicCircuit(double sampleRate);

    void addComponent(std::shared_ptr<Component> c);
    void addPotentiometer(std::shared_ptr<Potentiometer> p);
    void setSourceImpedance(double impedance);

    double process(double inputVoltage, int probeNode);

    // Get a potentiometer by index for real-time adjustment
    std::shared_ptr<Potentiometer> getPotentiometer(int index);
};

} // namespace MKAudio
