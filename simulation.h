#pragma once

#include <memory>
#include <array>

namespace MKAudio
{

// ==========================================
// 1. Optimized Component Classes
// ==========================================
class Component {
protected:
    int nodeA, nodeB;

public:
    Component(int n1, int n2) : nodeA(n1), nodeB(n2) {}
    virtual ~Component() = default;

    // 1. Preprocess Step: Contributes to the static Y Matrix
    virtual double getConductance(double dt) const = 0;

    // 2. Real-time Step: Contributes to the dynamic J Vector
    virtual double getCurrentSource(double dt) const = 0;

    // 3. Post-Process Step: Updates internal memory
    virtual void updateState(double vA, double vB, double dt) = 0;

    int getNodeA() const { return nodeA; }
    int getNodeB() const { return nodeB; }
};

class Resistor : public Component {
    double conductance;
public:
    Resistor(int n1, int n2, double r) : Component(n1, n2) { conductance = 1.0 / r; }
    
    double getConductance(double dt) const override { return conductance; }
    double getCurrentSource(double dt) const override { return 0.0; } // Resistive only
    void updateState(double vA, double vB, double dt) override {}
};

class Capacitor : public Component {
    double capacitance;
    double prevVoltage; 

public:
    Capacitor(int n1, int n2, double c) : Component(n1, n2), capacitance(c), prevVoltage(0.0) {}

    // G = C / dt
    double getConductance(double dt) const override { return capacitance / dt; }

    // J_eq = G * V_last
    double getCurrentSource(double dt) const override { return (capacitance / dt) * prevVoltage; }

    void updateState(double vA, double vB, double dt) override { prevVoltage = vA - vB; }
};

class Inductor : public Component {
    double inductance;
    double prevCurrent; 

public:
    Inductor(int n1, int n2, double l) : Component(n1, n2), inductance(l), prevCurrent(0.0) {}

    // G = dt / L
    double getConductance(double dt) const override { return dt / inductance; }

    // J_eq = -I_last
    double getCurrentSource(double dt) const override { return -prevCurrent; }

    void updateState(double vA, double vB, double dt) override {
        // Trapezoidal or Backward Euler update for current
        double voltage = vA - vB;
        prevCurrent += (voltage * dt) / inductance;
    }
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
    Circuit(double sampleRate) : currentDevices(0) { dt = 1.0 / sampleRate; }

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

} // namespace MKAudio