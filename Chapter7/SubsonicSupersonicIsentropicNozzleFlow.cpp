#include "SubsonicSupersonicIsentropicNozzleFlow.h"
#include <cmath>
#include <algorithm>

using namespace std;

SubsonicSupersonic::SubsonicSupersonic(double deltax, std::size_t N)
    : m_DeltaX(deltax), m_InvDeltaX(1.0 / deltax), m_N(N),
    m_X(N + 1),
    m_A(N + 1),
    m_LnA(N + 1),
    m_LnAGradForward(N + 1),
    m_LnAGradBackward(N + 1),
    m_DeltaTVec(N + 1),
    m_Rho(N + 1),
    m_V(N + 1),
    m_T(N + 1),
    m_RhoEstimated(N + 1),
    m_VEstimated(N + 1),
    m_TEstimated(N + 1),
    m_RhoGradEstimated(N + 1),
    m_VGradEstimated(N + 1),
    m_TGradEstimated(N + 1),
    m_RhoGradCorrected(N + 1),
    m_VGradCorrected(N + 1),
    m_TGradCorrected(N + 1)
{
    for (size_t i = 0; i < N + 1; ++i)
    {
        m_X[i] = i * m_DeltaX;
        m_A[i] = GetA(m_X[i]);
        m_LnA[i] = log(m_A[i]);
        m_Rho[i] = 1 - 0.3146 * m_X[i];
        m_T[i] = 1 - 0.2314 * m_X[i];
        m_V[i] = (0.1 + 1.09 * m_X[i]) * sqrt(m_T[i]);
    }
    for (size_t i = 1; i < N; ++i)
    {
        m_LnAGradForward[i] = m_InvDeltaX * (m_LnA[i + 1] - m_LnA[i]);
        m_LnAGradBackward[i] = m_InvDeltaX * (m_LnA[i] - m_LnA[i - 1]);
    }
}

void SubsonicSupersonic::Calculate(double gamma, double C, std::size_t iterate_times, std::ostream &out)
{
    double invgamma = 1.0 / gamma;
    double gamma_1 = gamma - 1.0;
    m_Rho[0] = m_RhoEstimated[0] = 1.0;
    m_T[0] = m_TEstimated[0] = 1.0;
    for (size_t cnt = 0; cnt < iterate_times; ++cnt)
    {
        for (size_t i = 1; i < m_N; ++i)
        {
            m_DeltaTVec[i] = C * m_DeltaX / (m_V[i] + sqrt(m_T[i]));
        }
        m_DeltaT = *min_element(m_DeltaTVec.begin() + 1, m_DeltaTVec.end() - 1);
        double halfDeltaT = 0.5 * m_DeltaT;

        for (size_t i = 1; i < m_N; ++i)
        {
            double Rhotmp = m_InvDeltaX * (m_Rho[i + 1] - m_Rho[i]);
            double Vtmp = m_InvDeltaX * (m_V[i + 1] - m_V[i]);
            double Ttmp = m_InvDeltaX * (m_T[i + 1] - m_T[i]);

            m_RhoGradEstimated[i] = -m_Rho[i] * Vtmp - m_V[i] * Rhotmp - m_Rho[i] * m_V[i] * m_LnAGradForward[i];
            m_VGradEstimated[i] = -m_V[i] * Vtmp - invgamma * (Ttmp + Rhotmp * m_T[i] / m_Rho[i]);
            m_TGradEstimated[i] = -m_V[i] * Ttmp - gamma_1 * m_T[i] * (Vtmp + m_V[i] * m_LnAGradForward[i]);
  
            m_RhoEstimated[i] = m_Rho[i] + m_RhoGradEstimated[i] * m_DeltaT;
            m_VEstimated[i] = m_V[i] + m_VGradEstimated[i] * m_DeltaT;
            m_TEstimated[i] = m_T[i] + m_TGradEstimated[i] * m_DeltaT;
        }
        m_VEstimated[0] = 2 * m_VEstimated[1] - m_VEstimated[2];
        m_RhoEstimated[m_N] = 2 * m_RhoEstimated[m_N - 1] - m_RhoEstimated[m_N - 2];
        m_VEstimated[m_N] = 2 * m_VEstimated[m_N - 1] - m_VEstimated[m_N - 2];
        m_TEstimated[m_N] = 2 * m_TEstimated[m_N - 1] - m_TEstimated[m_N - 2];

        for (size_t i = 1; i < m_N; ++i)
        {
            double Rhotmp = m_InvDeltaX * (m_RhoEstimated[i] - m_RhoEstimated[i - 1]);
            double Vtmp = m_InvDeltaX * (m_VEstimated[i] - m_VEstimated[i - 1]);
            double Ttmp = m_InvDeltaX * (m_TEstimated[i] - m_TEstimated[i - 1]);

            m_RhoGradCorrected[i] = -m_RhoEstimated[i] * Vtmp - m_VEstimated[i] * Rhotmp - m_RhoEstimated[i] * m_VEstimated[i] * m_LnAGradBackward[i];
            m_VGradCorrected[i] = -m_VEstimated[i] * Vtmp - invgamma * (Ttmp + Rhotmp * m_TEstimated[i] / m_RhoEstimated[i]);
            m_TGradCorrected[i] = -m_VEstimated[i] * Ttmp - gamma_1 * m_TEstimated[i] * (Vtmp + m_VEstimated[i] * m_LnAGradBackward[i]);

            m_Rho[i] += halfDeltaT * (m_RhoGradEstimated[i] + m_RhoGradCorrected[i]);
            m_V[i] += halfDeltaT * (m_VGradEstimated[i] + m_VGradCorrected[i]);
            m_T[i] += halfDeltaT * (m_TGradEstimated[i] + m_TGradCorrected[i]);
        }
        m_V[0] = 2 * m_V[1] - m_V[2];
        m_Rho[m_N] = 2 * m_Rho[m_N - 1] - m_Rho[m_N - 2];
        m_V[m_N] = 2 * m_V[m_N - 1] - m_V[m_N - 2];
        m_T[m_N] = 2 * m_T[m_N - 1] - m_T[m_N - 2];
    }
    out << "i" << "," << "x" << "," << "Rho" << "," << "V" << "," << "T" << "\n";
    for (size_t i = 0; i < m_N + 1; ++i)
    {
        out << i + 1 << "," << m_X[i] << "," << m_Rho[i] << "," << m_V[i] << "," << m_T[i] << "\n";
    }
}
