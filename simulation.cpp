#include <cmath>
#include <algorithm>

#include "simulation.h"

namespace MKAudio
{

template <std::size_t NumDevices, std::size_t NumNodes>
void Circuit<NumDevices, NumNodes>::solveLinearSystem(int N)
{
    // Copy static Y to work Y because elimination destroys the matrix
    std::copy(Y_static.begin(), Y_static.begin() + N * NumNodes, Y_work.begin());
    // Note: We copy N rows. Stride is still NumNodes.

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
            // Optimization: Start inner loop from i to avoid zeroing out known columns
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
void Circuit<NumDevices, NumNodes>::addComponent(std::shared_ptr<Component> c) { if(currentDevices < NumDevices) components[++currentDevices] = c; }

template <std::size_t NumDevices, std::size_t NumNodes>
void Circuit<NumDevices, NumNodes>::preprocess(double impedence)
{
    // Clear Matrix
    Y_static.fill(0.0);

    for (std::size_t i = 0; i < NumDevices; ++i) {
        auto& comp = components[i];
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
    // This ensures the matrix is non-singular even if the user just attaches a capacitor
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
    for (std::size_t i = 0; i < NumDevices; ++i) {
        auto& comp = components[i];
        double Is = comp->getCurrentSource(dt); // L/C history
        if (Is == 0.0) continue;

        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();

        // Current enters node A, leaves node B (or vice versa depending on definition)
        // If component is defined A->B, and Is is positive:
        // It injects current OUT of A and INTO B?
        // Standard MNA: Right Hand Side vector J accumulates KNOWN currents entering the node.
        // Our Companion models usually calculate Is parallel to G.
        // Check signs: For Cap, Is = (C/dt)*V_old. This acts as a source "pushing" back.

        if (n1 > 0) J[n1 - 1] -= Is;
        if (n2 > 0) J[n2 - 1] += Is;
    }

    // 4. Solve for Voltages
    solveLinearSystem(NumNodes);

    // 5. Update Component States
    for (std::size_t i = 0; i < NumDevices; ++i) {
        auto& comp = components[i];
        int n1 = comp->getNodeA();
        int n2 = comp->getNodeB();

        double v1 = (n1 == 0) ? 0.0 : nodes[n1 - 1];
        double v2 = (n2 == 0) ? 0.0 : nodes[n2 - 1];

        comp->updateState(v1, v2, dt);
    }

    if (probeNode <= 0 || probeNode > NumNodes) return 0.0;
    return nodes[probeNode - 1];
}

} // namespace MKAudio