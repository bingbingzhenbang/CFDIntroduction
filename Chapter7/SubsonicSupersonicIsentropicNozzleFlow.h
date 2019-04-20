#pragma once

#include <vector>
#include <iostream>

class SubsonicSupersonic
{
private:
    std::size_t m_N;
    double m_Gamma = 1.4;
    double m_C = 0.5;
    double m_DeltaX;
    double m_InvDeltaX;
    double m_DeltaT;

    std::vector<double> m_X;
    std::vector<double> m_A;
    std::vector<double> m_LnA;
    std::vector<double> m_LnAGradForward;
    std::vector<double> m_LnAGradBackward;
    std::vector<double> m_DeltaTVec;

    std::vector<double> m_Rho;
    std::vector<double> m_V;
    std::vector<double> m_T;

    std::vector<double> m_RhoEstimated;
    std::vector<double> m_VEstimated;
    std::vector<double> m_TEstimated;

    std::vector<double> m_RhoGradEstimated;
    std::vector<double> m_VGradEstimated;
    std::vector<double> m_TGradEstimated;

    std::vector<double> m_RhoGradCorrected;
    std::vector<double> m_VGradCorrected;
    std::vector<double> m_TGradCorrected;

public:
    double GetA(double x)
    {
        double tmp = x - 1.5;
        return 1 + 2.2 * tmp * tmp;
    }
    SubsonicSupersonic(double deltax, std::size_t N);
    void Calculate(double gamma, double C, std::size_t iterate_times, std::ostream &out);
};
