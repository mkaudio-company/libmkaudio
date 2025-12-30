#include <cmath>
#include <algorithm>

#include "simulation.h"

namespace MKAudio
{

// ==========================================
// Component Base Class
// ==========================================
Component::Component(int n1, int n2) : nodeA(n1), nodeB(n2) {}

bool Component::isDynamic() const { return false; }

int Component::getNodeA() const { return nodeA; }

int Component::getNodeB() const { return nodeB; }

// ==========================================
// Resistor
// ==========================================
Resistor::Resistor(int n1, int n2, double r) : Component(n1, n2) {
    conductance = 1.0 / r;
}

double Resistor::getConductance(double dt) const { return conductance; }

double Resistor::getCurrentSource(double dt) const { return 0.0; }

void Resistor::updateState(double vA, double vB, double dt) {}

void Resistor::setResistance(double r) { conductance = 1.0 / std::max(r, 1.0); }

double Resistor::getResistance() const { return 1.0 / conductance; }

// ==========================================
// Capacitor
// ==========================================
Capacitor::Capacitor(int n1, int n2, double c)
    : Component(n1, n2), capacitance(c), prevVoltage(0.0) {}

double Capacitor::getConductance(double dt) const {
    return capacitance / dt;
}

double Capacitor::getCurrentSource(double dt) const {
    return (capacitance / dt) * prevVoltage;
}

void Capacitor::updateState(double vA, double vB, double dt) {
    prevVoltage = vA - vB;
}

// ==========================================
// Inductor
// ==========================================
Inductor::Inductor(int n1, int n2, double l)
    : Component(n1, n2), inductance(l), prevCurrent(0.0) {}

double Inductor::getConductance(double dt) const {
    return dt / inductance;
}

double Inductor::getCurrentSource(double dt) const {
    return -prevCurrent;
}

void Inductor::updateState(double vA, double vB, double dt) {
    double voltage = vA - vB;
    prevCurrent += (voltage * dt) / inductance;
}

// ==========================================
// Switch
// ==========================================
void Switch::updateConductance() {
    conductance = closed ? (1.0 / R_ON) : (1.0 / R_OFF);
}

Switch::Switch(int n1, int n2, bool initialState)
    : Component(n1, n2), closed(initialState)
{
    updateConductance();
}

bool Switch::isDynamic() const { return true; }

double Switch::getConductance(double dt) const { return conductance; }

double Switch::getCurrentSource(double dt) const { return 0.0; }

void Switch::updateState(double vA, double vB, double dt) {}

void Switch::setClosed(bool state) {
    closed = state;
    updateConductance();
}

void Switch::setOpen(bool state) {
    closed = !state;
    updateConductance();
}

void Switch::toggle() {
    closed = !closed;
    updateConductance();
}

bool Switch::isClosed() const { return closed; }

bool Switch::isOpen() const { return !closed; }

// ==========================================
// Potentiometer
// ==========================================
void Potentiometer::updateConductances() {
    double rUpper = totalResistance * (1.0 - position);
    double rLower = totalResistance * position;

    // Clamp to avoid zero resistance
    rUpper = std::max(rUpper, MIN_R);
    rLower = std::max(rLower, MIN_R);

    conductanceUpper = 1.0 / rUpper;
    conductanceLower = 1.0 / rLower;
}

Potentiometer::Potentiometer(int n1, int n2, int nW, double totalR, double initialPosition)
    : Component(n1, n2), nodeWiper(nW), totalResistance(totalR), position(initialPosition)
{
    updateConductances();
}

bool Potentiometer::isDynamic() const { return true; }

double Potentiometer::getConductance(double dt) const {
    return conductanceUpper + conductanceLower;
}

double Potentiometer::getCurrentSource(double dt) const { return 0.0; }

void Potentiometer::updateState(double vA, double vB, double dt) {}

void Potentiometer::setPosition(double pos) {
    position = std::clamp(pos, 0.0, 1.0);
    updateConductances();
}

void Potentiometer::setPositionLog(double pos) {
    pos = std::clamp(pos, 0.0, 1.0);
    // Audio taper approximation: use logarithmic curve
    // At 50% rotation, output is ~10% of max
    double logPos = (std::exp(pos * 3.0) - 1.0) / (std::exp(3.0) - 1.0);
    position = logPos;
    updateConductances();
}

void Potentiometer::setPositionRevLog(double pos) {
    pos = std::clamp(pos, 0.0, 1.0);
    double revLogPos = std::log(1.0 + pos * (std::exp(3.0) - 1.0)) / 3.0;
    position = revLogPos;
    updateConductances();
}

double Potentiometer::getPosition() const { return position; }

int Potentiometer::getWiperNode() const { return nodeWiper; }

double Potentiometer::getConductanceUpper() const { return conductanceUpper; }

double Potentiometer::getConductanceLower() const { return conductanceLower; }

double Potentiometer::getTotalResistance() const { return totalResistance; }

// ==========================================
// Circuit Template Implementation
// ==========================================
template <std::size_t NumDevices, std::size_t NumNodes>
Circuit<NumDevices, NumNodes>::Circuit(double sampleRate) : currentDevices(0) {
    dt = 1.0 / sampleRate;
}

template <std::size_t NumDevices, std::size_t NumNodes>
void Circuit<NumDevices, NumNodes>::solveLinearSystem(int N)
{
    // Copy static Y to work Y because elimination destroys the matrix
    std::copy(Y_static.begin(), Y_static.begin() + N * NumNodes, Y_work.begin());

    for (int i = 0; i < N; ++i) {
        // Pivot
        int pivot = i;
        double maxVal = std::abs(Y_work[i * NumNodes + i]);

        for (int k = i + 1; k < N; ++k) {
            double val = std::abs(Y_work[k * NumNodes + i]);
            if (val > maxVal) {
                maxVal = val;
                pivot = k;
            }
        }

        // Swap rows in Matrix (Y) and Vector (J)
        if (pivot != i) {
            for (int col = i; col < N; ++col) {
                std::swap(Y_work[i * NumNodes + col], Y_work[pivot * NumNodes + col]);
            }
            std::swap(J[i], J[pivot]);
        }

        // Eliminate
        double pivotVal = Y_work[i * NumNodes + i];
        if (std::abs(pivotVal) < 1e-9) continue; // Singularity check

        for (int k = i + 1; k < N; ++k) {
            double factor = Y_work[k * NumNodes + i] / pivotVal;
            for (int j = i; j < N; ++j) {
                Y_work[k * NumNodes + j] -= factor * Y_work[i * NumNodes + j];
            }
            J[k] -= factor * J[i];
        }
    }

    // Back Substitution
    for (int i = N - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < N; ++j) {
            sum += Y_work[i * NumNodes + j] * nodes[j];
        }
        nodes[i] = (J[i] - sum) / Y_work[i * NumNodes + i];
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void Circuit<NumDevices, NumNodes>::addComponent(std::shared_ptr<Component> c) {
    if (currentDevices < static_cast<int>(NumDevices)) {
        components[currentDevices++] = c;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void Circuit<NumDevices, NumNodes>::preprocess(double impedence)
{
    // Clear Matrix
    Y_static.fill(0.0);

    for (int i = 0; i < currentDevices; ++i) {
        auto& comp = components[i];
        if (!comp) continue;

        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();
        double G = comp->getConductance(dt);

        // Stamp Y Matrix (0-indexed logic: Node 1 is index 0)
        if (n1 > 0) Y_static[(n1 - 1) * NumNodes + (n1 - 1)] += G;
        if (n2 > 0) Y_static[(n2 - 1) * NumNodes + (n2 - 1)] += G;

        if (n1 > 0 && n2 > 0) {
            Y_static[(n1 - 1) * NumNodes + (n2 - 1)] -= G;
            Y_static[(n2 - 1) * NumNodes + (n1 - 1)] -= G;
        }
    }

    // Add Source Resistance for Node 1 (Input)
    double G_source = impedence;
    if (NumNodes >= 1) {
        Y_static[0] += G_source;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
double Circuit<NumDevices, NumNodes>::process(double inputVoltage, int probeNode)
{
    // 1. Reset J Vector (Currents)
    std::fill(J.begin(), J.begin() + NumNodes, 0.0);

    // 2. Add Input Source (Norton Equivalent at Node 1)
    double G_source = 1.0 / 0.1;
    J[0] += inputVoltage * G_source;

    // 3. Accumulate Dynamic Currents from Components (Memory)
    for (int i = 0; i < currentDevices; ++i) {
        auto& comp = components[i];
        if (!comp) continue;

        double Is = comp->getCurrentSource(dt);
        if (Is == 0.0) continue;

        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();

        if (n1 > 0) J[n1 - 1] -= Is;
        if (n2 > 0) J[n2 - 1] += Is;
    }

    // 4. Solve for Voltages
    solveLinearSystem(NumNodes);

    // 5. Update Component States
    for (int i = 0; i < currentDevices; ++i) {
        auto& comp = components[i];
        if (!comp) continue;

        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();

        double v1 = (n1 == 0) ? 0.0 : nodes[n1 - 1];
        double v2 = (n2 == 0) ? 0.0 : nodes[n2 - 1];

        comp->updateState(v1, v2, dt);
    }

    if (probeNode <= 0 || probeNode > static_cast<int>(NumNodes)) return 0.0;
    return nodes[probeNode - 1];
}

// ==========================================
// DynamicCircuit Template Implementation
// ==========================================
template <std::size_t NumDevices, std::size_t NumNodes>
DynamicCircuit<NumDevices, NumNodes>::DynamicCircuit(double sampleRate) : dt(1.0 / sampleRate) {}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::solveLinearSystem(int N) {
    for (int i = 0; i < N; ++i) {
        int pivot = i;
        double maxVal = std::abs(Y_work[i * NumNodes + i]);

        for (int k = i + 1; k < N; ++k) {
            double val = std::abs(Y_work[k * NumNodes + i]);
            if (val > maxVal) {
                maxVal = val;
                pivot = k;
            }
        }

        if (pivot != i) {
            for (int col = i; col < N; ++col) {
                std::swap(Y_work[i * NumNodes + col], Y_work[pivot * NumNodes + col]);
            }
            std::swap(J[i], J[pivot]);
        }

        double pivotVal = Y_work[i * NumNodes + i];
        if (std::abs(pivotVal) < 1e-12) continue;

        for (int k = i + 1; k < N; ++k) {
            double factor = Y_work[k * NumNodes + i] / pivotVal;
            for (int j = i; j < N; ++j) {
                Y_work[k * NumNodes + j] -= factor * Y_work[i * NumNodes + j];
            }
            J[k] -= factor * J[i];
        }
    }

    for (int i = N - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < N; ++j) {
            sum += Y_work[i * NumNodes + j] * nodes[j];
        }
        double denom = Y_work[i * NumNodes + i];
        nodes[i] = (std::abs(denom) > 1e-12) ? (J[i] - sum) / denom : 0.0;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::stampComponent(std::shared_ptr<Component>& comp) {
    int n1 = comp->getNodeA();
    int n2 = comp->getNodeB();
    double G = comp->getConductance(dt);

    if (n1 > 0) Y_work[(n1 - 1) * NumNodes + (n1 - 1)] += G;
    if (n2 > 0) Y_work[(n2 - 1) * NumNodes + (n2 - 1)] += G;

    if (n1 > 0 && n2 > 0) {
        Y_work[(n1 - 1) * NumNodes + (n2 - 1)] -= G;
        Y_work[(n2 - 1) * NumNodes + (n1 - 1)] -= G;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::stampPotentiometer(std::shared_ptr<Potentiometer>& pot) {
    int nA = pot->getNodeA();
    int nB = pot->getNodeB();
    int nW = pot->getWiperNode();

    double Gu = pot->getConductanceUpper(); // A to Wiper
    double Gl = pot->getConductanceLower(); // Wiper to B

    // Stamp upper resistor (nA to nW)
    if (nA > 0) Y_work[(nA - 1) * NumNodes + (nA - 1)] += Gu;
    if (nW > 0) Y_work[(nW - 1) * NumNodes + (nW - 1)] += Gu;
    if (nA > 0 && nW > 0) {
        Y_work[(nA - 1) * NumNodes + (nW - 1)] -= Gu;
        Y_work[(nW - 1) * NumNodes + (nA - 1)] -= Gu;
    }

    // Stamp lower resistor (nW to nB)
    if (nW > 0) Y_work[(nW - 1) * NumNodes + (nW - 1)] += Gl;
    if (nB > 0) Y_work[(nB - 1) * NumNodes + (nB - 1)] += Gl;
    if (nW > 0 && nB > 0) {
        Y_work[(nW - 1) * NumNodes + (nB - 1)] -= Gl;
        Y_work[(nB - 1) * NumNodes + (nW - 1)] -= Gl;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::addComponent(std::shared_ptr<Component> c) {
    if (currentDevices < static_cast<int>(NumDevices)) {
        components[currentDevices++] = c;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::addPotentiometer(std::shared_ptr<Potentiometer> p) {
    if (currentPots < static_cast<int>(NumDevices)) {
        potentiometers[currentPots++] = p;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::addSwitch(std::shared_ptr<Switch> s) {
    if (currentSwitches < static_cast<int>(NumDevices)) {
        switches[currentSwitches++] = s;
    }
}

template <std::size_t NumDevices, std::size_t NumNodes>
void DynamicCircuit<NumDevices, NumNodes>::setSourceImpedance(double impedance) {
    sourceImpedance = impedance;
}

template <std::size_t NumDevices, std::size_t NumNodes>
double DynamicCircuit<NumDevices, NumNodes>::process(double inputVoltage, int probeNode) {
    // 1. Clear Y matrix and J vector
    Y_work.fill(0.0);
    std::fill(J.begin(), J.begin() + NumNodes, 0.0);

    // 2. Stamp all static components
    for (int i = 0; i < currentDevices; ++i) {
        if (components[i]) {
            stampComponent(components[i]);
        }
    }

    // 3. Stamp all potentiometers (dynamic)
    for (int i = 0; i < currentPots; ++i) {
        if (potentiometers[i]) {
            stampPotentiometer(potentiometers[i]);
        }
    }

    // 4. Stamp all switches (dynamic) - switches are simple 2-terminal components
    for (int i = 0; i < currentSwitches; ++i) {
        if (switches[i]) {
            // Switches can be stamped like regular components since they're 2-terminal
            std::shared_ptr<Component> switchComp = switches[i];
            stampComponent(switchComp);
        }
    }

    // 5. Add source impedance at node 1
    double G_source = 1.0 / sourceImpedance;
    if (NumNodes >= 1) {
        Y_work[0] += G_source;
    }

    // 6. Add input current source (Norton equivalent)
    J[0] += inputVoltage * G_source;

    // 7. Accumulate dynamic currents from reactive components
    for (int i = 0; i < currentDevices; ++i) {
        auto& comp = components[i];
        if (!comp) continue;

        double Is = comp->getCurrentSource(dt);
        if (Is == 0.0) continue;

        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();

        if (n1 > 0) J[n1 - 1] -= Is;
        if (n2 > 0) J[n2 - 1] += Is;
    }

    // 8. Solve linear system
    solveLinearSystem(NumNodes);

    // 9. Update component states
    for (int i = 0; i < currentDevices; ++i) {
        auto& comp = components[i];
        if (!comp) continue;

        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();

        double v1 = (n1 == 0) ? 0.0 : nodes[n1 - 1];
        double v2 = (n2 == 0) ? 0.0 : nodes[n2 - 1];

        comp->updateState(v1, v2, dt);
    }

    if (probeNode <= 0 || probeNode > static_cast<int>(NumNodes)) return 0.0;
    return nodes[probeNode - 1];
}

template <std::size_t NumDevices, std::size_t NumNodes>
std::shared_ptr<Potentiometer> DynamicCircuit<NumDevices, NumNodes>::getPotentiometer(int index) {
    if (index >= 0 && index < currentPots) {
        return potentiometers[index];
    }
    return nullptr;
}

template <std::size_t NumDevices, std::size_t NumNodes>
std::shared_ptr<Switch> DynamicCircuit<NumDevices, NumNodes>::getSwitch(int index) {
    if (index >= 0 && index < currentSwitches) {
        return switches[index];
    }
    return nullptr;
}

// ==========================================
// Explicit Template Instantiations
// ==========================================
// Circuit types used in the project
template class Circuit<3, 2>;   // Power amp circuit

// DynamicCircuit types used in the project
template class DynamicCircuit<10, 5>;  // Fender tone stack

} // namespace MKAudio
