/******************************************************************************
 *   Copyright (C) 2025       						      *
 *               - Aldric Revel, MSU, USA                                     *
 *               - Jose Luis Rodriguez-Sánchez, University of Coruña, Spain   *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 ******************************************************************************/

#include "TApplication.h"
#include "mass.h"
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStreamerInfo.h>
#include <TStyle.h>
#include <TTree.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <thread>
#include <vector>

void showProgressBar(int i, int total, int width = 50)
{
    float ratio = i / (total - 1.);
    int c = int(ratio * width);

    std::cout << "Processed: \033[32m" << i + 1 << "\033[0m / \033[31m" << total << "\033[0m [";

    for (int j = 0; j < width; j++)
    {
        if (j < c)
            std::cout << "\033[32m█\033[0m";
        else
            std::cout << "-";
    }

    std::cout << "] \033[33m" << int(ratio * 100) << "%\033[0m" << '\r' << std::flush;
}

double VMOD(double* vec)
{
    double res = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    return res;
}

double VoV(double* V, double* W)
{ // scalar product
    double res = V[0] * W[0] + V[1] * W[1] + V[2] * W[2];
    return res;
}

double Cnn_0(double X, double r0)
{ // For nn-FSI in direct decay

    double hc, f0, d0;
    hc = 197.3; //[MeVfm]
    f0 = 18.5;  // nn scatt. length [fm]
    d0 = 2.8;   // effective range  [fm]
    double a, b, Ref, Imf, f_2, B0, Bi, pi, F1, F2;
    double res = -40000;

    if (r0 == 0.0)
    {
        res = 1.;
    }
    else if (r0 < 1.1)
    { // modified on 11th december 2019 (initial value: r0<1.5)
        cout << "r0 too small" << endl;
    }
    else
    {
        pi = acos(-1.);
        a = 1. / f0 + d0 * (X / hc) * (X / hc) / 2.;
        b = -X / hc; // f = 1/(a+ib)
        Ref = a / (a * a + b * b);
        Imf = -b / (a * a + b * b);
        f_2 = Ref * Ref + Imf * Imf;
        B0 = -0.5 * exp(-4. * (r0 * X / hc) * (r0 * X / hc)); // t0 = 0

        if (X == 0.)
        {
            F1 = 1.;
            F2 = 0.;
        }
        else
        {
            F1 = (TMath::Erf(2. * r0 * X / hc) * sqrt(pi) / 2.) * exp(-4. * (r0 * X / hc) * (r0 * X / hc)) /
                 (2. * r0 * X / hc);
            F2 = (1. - exp(-4. * (r0 * X / hc) * (r0 * X / hc))) / (2. * r0 * X / hc);
        }
        Bi = 0.25 * f_2 / (r0 * r0) * (1. - d0 / (2. * sqrt(pi) * r0)) + Ref * F1 / (sqrt(pi) * r0) -
             Imf * F2 / (2. * r0);
        res = 1. + B0 + Bi;
    }
    return res;
}

double intBi(double Y,
             double r0,
             double t0,
             double rT,
             double G,
             double Mv,
             double rho,
             double f_2,
             double k,
             double Imf,
             double qT,
             double q0,
             double Ref)
{ // (X,Y) = (rT,rL)

    double hc, rL, r;
    hc = 197.3; // hc [MeVfm]

    rL = Y;
    r = TMath::Sqrt(rT * rT + rL * rL);
    // A  = (r0*t0)/(G*TMath::Sqrt(t0*t0+(Mv*r0)*(Mv*r0)));

    double res =
        rT * TMath::Exp(-(rT / (2. * r0)) * (rT / (2. * r0)) - (rL / (2. * G * rho)) * (rL / (2. * G * rho))) *
        (f_2 / (2. * r * r) + (Ref * cos(k * r / hc) - Imf * TMath::Sin(k * r / hc)) / r *
                                  TMath::BesselJ0(.5 * qT * rT / hc) * TMath::Cos(q0 * rL / (2. * G * Mv * hc)));

    return res;
}

double Integrate_BI(double r0,
                    double t0,
                    double G,
                    double Mv,
                    double rho,
                    double f_2,
                    double k,
                    double Imf,
                    double qT,
                    double q0,
                    double Ref)
{
    double intXY, step, x1, xn, x2, intX, rT;
    int i, j, N;
    step = 1.; // 0.05; //size of the integration step
    x1 = 0.;   // start of the integration
    x2 = 20.;  // 5. //end of the integration

    intXY = 0.;
    N = int((x2 - x1) / step);

    for (i = 1; i <= N; i++)
    {
        rT = i * step - step / 2.;
        intX = 0.;

        for (j = 1; j <= N; j++)
        {
            xn = x1 + j * step - step / 2.;
            intX = intX + intBi(xn, r0, t0, rT, G, Mv, rho, f_2, k, Imf, qT, q0, Ref) * step;
        }
        intXY = intXY + intX * step;
    }
    return intXY;
}

double Cnn_sequential(double* p1, double* p2, double r0, double t0)
{ // eq (21)

    int i;
    double hc, f0, d0;
    hc = 197.3; //[MeVfm]
    f0 = 18.5;  // nn scatt. length [fm]
    d0 = 2.8;   // effective range  [fm]
    double a, b, Ref, Imf, f_2, B0, Bi, pi, integralBi, intWd, q, q0, k, Mv, G, qT, rho, u, C, mu, nu, tab;
    double res = -40000;
    double v[3];

    if (r0 == 0.)
    {
        res = 1.;
    }
    else if (r0 < 1.1)
    { //  modified on 11th december 2019 (initial value: r0<1.5)
        std::cout << "r0 too small" << std::endl;
    }
    else
    {
        pi = acos(-1.);
        q = TMath::Sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) +
                        (p1[2] - p2[2]) * (p1[2] - p2[2]));
        q0 = TMath::Abs(p1[3] - p2[3]);
        k = TMath::Sqrt(TMath::Abs(q * q - q0 * q0)) / 2.;

        for (i = 0; i < 3; i++)
        {
            v[i] = (p1[i] + p2[i]) / (p1[3] + p2[3]);
        }

        Mv = TMath::Sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); //"v" in article
        G = 1. / TMath::Sqrt(TMath::Abs(1. - Mv * Mv));            // gamma
        qT = TMath::Sqrt(TMath::Abs(4. * k * k - (q0 / (G * Mv)) * (q0 / (G * Mv))));
        rho = TMath::Sqrt(r0 * r0 + (Mv * t0) * (Mv * t0));
        a = 1. / f0 + d0 * (k / hc) * (k / hc) / 2.;
        b = -k / hc;
        Ref = a / (a * a + b * b);
        Imf = -b / (a * a + b * b);
        f_2 = Ref * Ref + Imf * Imf;

        B0 = -0.5 * exp(-(r0 * q / hc) * (r0 * q / hc) - (t0 * q0 / hc) * (t0 * q0 / hc));
        // cout << "Sequ =" << pi << "   " << a << "   " << b << "  " << Ref << " "
        // << Imf << "   " << f_2 << "   " << B0 << endl;
        integralBi = Integrate_BI(r0, t0, G, Mv, rho, f_2, k, Imf, qT, q0, Ref);

        Bi = integralBi / (2. * TMath::Sqrt(pi) * r0 * r0 * G * rho) -
             f_2 * d0 / (8. * TMath::Sqrt(pi) * G * rho * r0 * r0);

        res = 1. + B0 + Bi;
    }
    return res;
}

double Estar1n(double mn, double mf, double* Pnn1, double* Pnn3)
{ // Calculate Estar if only one neutron is
  // seen in the detector
    double res;
    res = TMath::Sqrt(mf * mf + mn * mn + 2 * Pnn1[3] * Pnn3[3] - 2 * VoV(Pnn1, Pnn3)) - mf - mn;
    return res;
}

double Estar(double mn, double mf, double* Pnn1, double* Pnn2, double* Pnn3)
{ // Calculate Estar if the 2 neutrons are seen as one
    double res;
    double Pnn[4];
    Pnn[0] = (Pnn1[0] + Pnn2[0]) / 2;
    Pnn[1] = (Pnn1[1] + Pnn2[1]) / 2;
    Pnn[2] = (Pnn1[2] + Pnn2[2]) / 2;
    Pnn[3] = Pnn1[3] + Pnn2[3] - mn;

    res = TMath::Sqrt(mf * mf + mn * mn + 2 * Pnn[3] * Pnn3[3] - 2 * VoV(Pnn, Pnn3)) - mf - mn;
    return res;
}

double Estar2n(double mn, double mf, double* Pnn1, double* Pnn2, double* Pnn3)
{
    double res;
    res = TMath::Sqrt(mf * mf + 2 * mn * mn + 2 * Pnn1[3] * Pnn3[3] - 2 * VoV(Pnn1, Pnn3) + 2 * Pnn2[3] * Pnn3[3] -
                      2 * VoV(Pnn2, Pnn3) + 2 * Pnn2[3] * Pnn1[3] - 2 * VoV(Pnn2, Pnn1)) -
          mf - 2 * mn;
    return res;
}

void ZBOOST(double* P, double beta)
{
    double gamma, Ppar, Etot;
    gamma = 1. / sqrt(1. - beta * beta);
    Ppar = P[2];
    Etot = P[3];
    P[2] = gamma * (Ppar + beta * Etot);
    P[3] = gamma * (Etot + beta * Ppar);
}

void BOOST(double* P, double* beta)
{
    // Variables
    int i;
    double mB, gamma, BoP, Ppar, Etot;
    double Pper[3];

    // BOOST
    mB = VMOD(beta);

    if (fabs(mB) > 1.e-6)
    {
        gamma = 1. / sqrt(1. - mB * mB);
        BoP = VoV(beta, P);
        Ppar = BoP / mB;
        Etot = P[3];

        for (i = 0; i < 3; i++)
        {
            Pper[i] = P[i] - Ppar * beta[i] / mB;
        }

        Ppar = gamma * (Ppar + (mB * Etot));
        P[3] = gamma * (Etot + BoP);

        for (i = 0; i < 3; i++)
        {
            P[i] = Pper[i] + Ppar * beta[i] / mB;
        }
    }
    else
    { /*cout << "no BOOST" << endl;*/
    }
}

double dNdp3(double a1, double a2, double a3, double E, double p3)
{
    double E3, res;
    E3 = sqrt(a3 * a3 + p3 * p3);
    res = (p3 * p3 / E3) *
          sqrt(fabs((E * E + a3 * a3 - 2. * E * E3 - (a1 + a2) * (a1 + a2)) *
                    (E * E + a3 * a3 - 2. * E * E3 - (a1 - a2) * (a1 - a2)))) /
          (E * E + a3 * a3 - 2. * E * E3);
    return res;
}

void PHASE_3(double a1, double a2, double a3, double Ed, double* w1, double* w2, double* w3, TRandom* Rand_angle)
{
    // Variables
    int Ndiv = 200;
    double ran[2];
    double beta[3];
    double E, p3M, M12, p3, p1, Pold, Fold, Pmax, Fmax, F3, step, Ec2n, Enn, the, phi;

    // 3rd particle + 2-body phase space
    E = a1 + a2 + a3 + Ed;
    p3M = sqrt((E * E - (a1 + a2 + a3) * (a1 + a2 + a3)) * (E * E - (a1 + a2 - a3) * (a1 + a2 - a3))) / (2. * E);
    Fold = 0.;
    step = p3M / Ndiv;
    Pmax = p3M - step;
    Fmax = dNdp3(a1, a2, a3, E, Pmax);
    while (Fmax > Fold)
    {
        Pold = Pmax;
        Fold = Fmax;
        Pmax = Pold - step;
        Fmax = dNdp3(a1, a2, a3, E, Pmax);
    }

    Fmax = dNdp3(a1, a2, a3, E, Pmax + step / 2.);

    do
    {
        ran[0] = Rand_angle->Uniform(0., 1.);
        ran[1] = Rand_angle->Uniform(0., 1.);
        p3 = p3M * ran[0];
        F3 = Fmax * ran[1];
    } while (F3 > dNdp3(a1, a2, a3, E, p3));

    ran[0] = Rand_angle->Uniform(0., 1.);
    ran[1] = Rand_angle->Uniform(0., 1.);
    the = acos(1. - 2. * ran[0]);
    phi = 2. * TMath::Pi() * ran[1];
    w3[0] = p3 * sin(the) * cos(phi);
    w3[1] = p3 * sin(the) * sin(phi);
    w3[2] = p3 * cos(the);
    w3[3] = sqrt(a3 * a3 + p3 * p3);

    // 2-body phase space in (12) cm
    M12 = sqrt((E - w3[3]) * (E - w3[3]) - p3 * p3);
    p1 = sqrt(fabs((M12 * M12 - (a1 + a2) * (a1 + a2)) * (M12 * M12 - (a1 - a2) * (a1 - a2)))) / (2. * M12);
    ran[0] = Rand_angle->Uniform(0., 1.);
    ran[1] = Rand_angle->Uniform(0., 1.);
    the = acos(1. - 2. * ran[0]);
    phi = 2. * TMath::Pi() * ran[1];
    w1[0] = p1 * sin(the) * cos(phi);
    w1[1] = p1 * sin(the) * sin(phi);
    w1[2] = p1 * cos(the);
    w1[3] = sqrt(a1 * a1 + p1 * p1);

    for (int i = 0; i < 3; i++)
    {
        w2[i] = -w1[i];
    }
    w2[3] = sqrt(a2 * a2 + p1 * p1);

    // Back to (123)CM
    for (int i = 0; i < 3; i++)
    {
        beta[i] = -w3[i] / (E - w3[3]);
    }
    BOOST(w2, beta);
    BOOST(w1, beta);
}

double R_2(double a1, double a2, double Ed, double Er)
{

    double m1, m2, Ex, E2, E, R_0;
    double res;

    m1 = a1;
    m2 = a2;
    Ex = Ed;
    E2 = Er;
    if (E2 >= Ex)
    {
        res = 0.;
    }
    else
    {
        E = m1 + m2 + Ex;
        R_0 = TMath::Sqrt((E * E - (m2 - m1) * (m2 - m1)) * (E * E - (m2 + m1) * (m2 + m1))) / (E * E);
        res =
            TMath::Sqrt((E * E - ((m2 + E2) - m1) * ((m2 + E2) - m1)) * (E * E - ((m2 + E2) + m1) * ((m2 + E2) + m1))) /
            (E * E) / R_0;
    }
    return res;
}

double PHASE2x2(double a1,
                double a2,
                double a3,
                double Ed,
                double* w1,
                double* w2,
                double* w3,
                TRandom* Rand_angle,
                double A,
                TF1* BW)
{ // sequential decay (give back resonance energy selected

    int i, l;
    double beta[3];
    double ran[2];
    double E23, the, phi;
    double m1, m2, Ex, M12, Pcm;
    double v1[4], v2[4];

    // Check first phase space......:
    l = 0;
    do
    {

        E23 = BW->GetRandom(0., Ed);
        l = l + 1;
        if (l > 50)
        {
            PHASE_3(a1, a2, a3, Ed, w1, w2, w3, Rand_angle);
            cout << "Too long to find proper Er" << endl;
            return 0.;
        }
        ran[0] = Rand_angle->Uniform(0., 1.);
    } while (ran[0] > R_2(a1, a2 + a3 + E23, Ed, E23) || E23 <= 0);

    // Go to DBLE......:
    m1 = a1;
    m2 = a2 + a3 + E23;
    Ex = Ed - E23;

    // 2-body phase-space (1+23)......:
    M12 = m1 + m2 + Ex;
    Pcm =
        TMath::Sqrt(TMath::Abs((M12 * M12 - (m1 + m2) * (m1 + m2)) * (M12 * M12 - (m1 - m2) * (m1 - m2)))) / (2. * M12);
    ran[0] = Rand_angle->Uniform(0., 1.);
    ran[1] = Rand_angle->Uniform(0., 1.);
    the = TMath::ACos(1. - 2. * ran[0]);
    phi = 2. * TMath::Pi() * ran[1];
    v1[0] = Pcm * sin(the) * cos(phi);
    v1[1] = Pcm * sin(the) * sin(phi);
    v1[2] = Pcm * cos(the);
    v1[3] = sqrt(m1 * m1 + Pcm * Pcm);
    for (i = 0; i < 3; i++)
    {
        v2[i] = -v1[i];
    }
    v2[3] = sqrt(m2 * m2 + Pcm * Pcm);
    for (i = 0; i < 4; i++)
    {
        w1[i] = v1[i];
    }

    //(23) 2-body decay......:
    for (i = 0; i < 3; i++)
    {
        beta[i] = v2[i] / v2[3];
    }
    m1 = a2;
    m2 = a3;
    Ex = E23;
    M12 = m1 + m2 + Ex;
    Pcm = sqrt(TMath::Abs((M12 * M12 - (m1 + m2) * (m1 + m2)) * (M12 * M12 - (m1 - m2) * (m1 - m2)))) / (2. * M12);
    ran[0] = Rand_angle->Uniform(0., 1.);
    ran[1] = Rand_angle->Uniform(0., 1.);
    the = TMath::ACos(1. - 2. * ran[0]);
    phi = 2. * TMath::Pi() * ran[1];
    v1[0] = Pcm * sin(the) * cos(phi);
    v1[1] = Pcm * sin(the) * sin(phi);
    v1[2] = Pcm * cos(the);
    v1[3] = sqrt(m1 * m1 + Pcm * Pcm);
    for (i = 0; i < 3; i++)
    {
        v2[i] = -v1[i];
    }
    v2[3] = sqrt(m2 * m2 + Pcm * Pcm);
    for (i = 0; i < 4; i++)
    {
        w2[i] = v1[i];
        w3[i] = v2[i];
    }
    BOOST(w2, beta);
    BOOST(w3, beta);

    return E23;
}

double DALITZ(double m_j, double m_i, double Ed, double* P_i, double* P_j)
{
    double res, M_ij;
    M_ij = m_i * m_i + m_j * m_j + 2 * P_i[3] * P_j[3] - 2 * VoV(P_i, P_j); // M_ij²
    res = (M_ij - (m_i + m_j) * (m_i + m_j)) / ((m_i + m_j + Ed) * (m_i + m_j + Ed) - (m_i + m_j) * (m_i + m_j));
    return res;
}

void PHASE_2(double a1, double a2, double Ed, double* w1, double* w2, TRandom* Rand_angle)
{
    int i;
    double the, phi, p1, Ex;
    double ran[2];

    Ex = a1 + a2 + Ed;
    p1 = sqrt(fabs((Ex * Ex - (a1 + a2) * (a1 + a2)) * (Ex * Ex - (a1 - a2) * (a1 - a2)))) / (2. * Ex);

    ran[0] = Rand_angle->Uniform(0., 1.);
    ran[1] = Rand_angle->Uniform(0., 1.);
    the = acos(1. - 2. * ran[0]);
    phi = 2. * TMath::Pi() * ran[1];

    w1[0] = p1 * sin(the) * cos(phi);
    w1[1] = p1 * sin(the) * sin(phi);
    w1[2] = p1 * cos(the);
    w1[3] = sqrt(a1 * a1 + p1 * p1);

    for (i = 0; i < 3; i++)
    {
        w2[i] = -w1[i];
    }
    w2[3] = sqrt(a2 * a2 + p1 * p1);
}

// A, Z of the fragment nuclei ex : 27Ne -> 26F -> 24F + 2n  SIMULATION(24,9)
void NeutronDecayGenerator(const TString Output_Name = "test",
                           const int Evt_number = 10,
                           const int A = 24,
                           const int Z = 9,
                           const double Ekin = 500., // MeV/u
                           const double E_BW = 2.5,
                           double W_BW = 0,
                           const int N = 1, // number of neutrons: 1n or 2n
                           const int decay_opt = 0)
{
    // Parameters of the Simulation  ---------------------------------------------------------------
    // E_BW is the resonance decay energy in MeV
    //
    // W_BW = 0 corresponds to a delta resonance with an average energy of E_BW in MeV
    // W_BW = -1 corresponds to an uniform distribution between 0 and E_BW in MeV
    // W_BW = -2 corresponds to a delta resonance with an average energy of 2 MeV
    // W_BW > 0 corresponds to a resonance with an average energy of E_BW in MeV and a sigma of W_BW

    std::cout << std::endl;
    std::cout << "Running for " << Evt_number << " events" << std::endl;
    std::cout << "Outgoing fragment A = " << A << " , Z = " << Z << std::endl;

    double Beta_sig = 0.; // Velocity dispersion for beta

    const double mf = masses_ame2020[A][Z]; // Mass of your fragment
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Outgoing fragment mass = " << mf <<" in MeV"<< std::endl;
    
    double r0_init = 1.5;                  // as defined in Lednicky paper
    const double Beta = TMath::Sqrt(1. - (mf / (Ekin * A + mf)) * (mf / (Ekin * A + mf)));
    int FSI_nn = 0;    // decay mode options
    int SEQ_Decay = 0; // decay mode options
    double E_DIN = 0.; // decay mode parameter
    double W_DIN = 0.; // decay mode parameter
    double Er_BW = 0.; // decay mode parameter
    double Wr_BW = 0.; // decay mode parameter
    const double fMax_energy_res = 12.;

    TString outfileName = Output_Name + (TString) ".out";
    ofstream outfile(outfileName, std::ofstream::out);

    /********PHASE SPACE DECAY************/
    if (decay_opt == 0)
    {
        std::cout << "Configuration for 1n decay" << std::endl;
    }
    else if (decay_opt == 1)
    {
        /********N-N CORRELATIONS*************/ //(Lednicky Paper)
        FSI_nn = 1;
        SEQ_Decay = 0;
        r0_init = 1.5; // should not be smaller than 1.5 according to paper and is
                       // link to rms distance between neutrons
    }
    else if (decay_opt == 2)
    {
        /*********SEQUENTIAL DECAY************/
        std::cout << "Configuration for sequential decay" << std::endl;
        FSI_nn = 1;
        SEQ_Decay = 1;
        r0_init = 1.5;
        Er_BW = 1.;
        Wr_BW = 0.1;
    }
    else if (decay_opt == 3)
    {
        /**********DINEUTRON DECAY************/
        std::cout << "Configuration for dineutron decay" << std::endl;
        FSI_nn = 1;
        SEQ_Decay = 2;
        E_DIN = 1.;
        W_DIN = 0.1;
    }
    else
    {
        std::cout << "This decay mode does not exist." << std::endl;
        return;
    }

    // Variables declaration
    const double amu = 931.494028;            // MeV/c2
    const double c = 299792458;               // m/s
    const double hc = 197.3269788;            // MeV.fm
    const double h_barC_MeV_fm = 197.3269788; // MeV.fm
    int Nb = Evt_number;
    double a1, a2, a3;
    double Pnn1[4], Pnn2[4], Pnn3[4]; // 4-vector with 2n
    double Pn1[4], Pn2[4];            // 4-vector with 1n
    double Pn[4], Pf[4];              // 4-vector with 1n
    double mn = 939.56542;            // MeV/c2

    double r0, t0, Gam_res;
    double Cnn_max; //= Cnn_0(0.,r0); //max value of the C_nn correlation in
                    // order to normalize it from 0 to 1
    double Er, E23, C_nn_seq;

    char* file_name_1n;
    char* file_name_2n;
    int f1;
    int f2;

    double b_beam, gamma_beam;

    // BRANCH VARIABLES
    double Ed;
    double Pn_X;
    double Pn_Y;
    double Pn_Z;
    double Pn_E;
    double Pf_X;
    double Pf_Y;
    double Pf_Z;
    double Pf_E;
    double Estar_1n;

    std::vector<double> m_nn;
    std::vector<double> m_fn1;
    std::vector<double> m_fn2;
    std::vector<double> E_star2n;

    std::vector<double> E_star2n_n1;
    std::vector<double> E_star2n_n2;
    std::vector<double> theta_n_n;
    std::vector<double> theta_fn1_n; // angle between p_fn1 and pn2
    std::vector<double> theta_fn2_n; // angle between p_fn2 and pn1
    std::vector<double> theta_nn1_f; // angle between pn1n2 and pf
    std::vector<double> theta_nn2_f; // angle between pn2n1 and pf
    std::vector<double> theta_kY_1;  // angle defined in 6Be paper k2->f k1->n1
    std::vector<double> theta_kT_1;  // angle defined in 6Be paper with k2->n2 k1->n1
    std::vector<double> Ex_Y_1;      // energy defined in 6Be paper k2->f k1->n1
    std::vector<double> Ex_T_1;      // energy defined in 6Be paper with k2->n2 k1->n1
    std::vector<double> theta_kY_2;  // angle defined in 6Be paper k2->f k1->n1
    std::vector<double> theta_kT_2;  // angle defined in 6Be paper with k2->n2 k1->n1
    std::vector<double> Ex_Y_2;      // energy defined in 6Be paper k2->f k1->n1
    std::vector<double> Ex_T_2;      // energy defined in 6Be paper with k2->n2 k1->n1
    double q;                        // argument to calculate before to use give it to C_nn

    double gamma_n[2], beta_n[2]; // gamma and beta of the two neutrons

    auto Rand_angle = std::make_unique<TRandom3>(0);
    auto Rand = std::make_unique<TRandom3>(0);
    auto Rand_0 = std::make_unique<TRandom3>(0);
    auto Dist = std::make_unique<TRandom3>(0);
    auto Dist1 = std::make_unique<TRandom3>(0);
    auto Dist2 = std::make_unique<TRandom3>(0);
    auto Rand_X1 = std::make_unique<TRandom3>(0);
    auto Rand_Y1 = std::make_unique<TRandom3>(0);
    auto Rand_X2 = std::make_unique<TRandom3>(0);
    auto Rand_Y2 = std::make_unique<TRandom3>(0);
    auto Rand_Z1 = std::make_unique<TRandom3>(0);
    auto Rand_Z2 = std::make_unique<TRandom3>(0);
    auto Time1 = std::make_unique<TRandom3>(0);
    auto Time2 = std::make_unique<TRandom3>(0);
    auto Rand_1n = std::make_unique<TRandom3>(0);
    auto Rand_2n = std::make_unique<TRandom3>(0);

    // angle_n_n in center of mass calculation
    double P_n1[4], P_n2[4], P_f[4], P_tot[4]; // 4-vector
    double P_n1n2[4], P_n2n1[4], P_fn1[4], P_fn2[4];
    double kx_Y_1[4], kx_T_1[4], ky_Y_1[4], ky_T_1[4];
    double kx_Y_2[4], kx_T_2[4], ky_Y_2[4], ky_T_2[4];
    double b_boost[3];

    // Histogram to visualize the nn-FSI correlation function
    // Histogram to visualize the nn-FSI correlation function
    auto h1_C_nn = std::make_unique<TH1D>("h1_C_nn", "h1_C_nn", 200, 0, 10);

    auto h1_C_nn_seq = std::make_unique<TH1D>("h1_C_nn_seq", "h1_C_nn_seq", 200, 0, 10); // from sequential decay

    auto h1_E_res = std::make_unique<TH1D>("h1_E_res", "h1_E_res", 200, 0, 15); // to check the intermediate resonance
                                                                                // energy in case of sequential decay

    // For simulation with Breit Wigner in Input
    double Xo_r = E_BW; // Centroid BW
    double Wo_r = W_BW; // GAMMA BW
    // Breit Wigner parameters
    double mu_r = 931.5 * A / (A + 1);
    double Ri_r = 1.4 * (pow(A, 1. / 3.) + pow(1, 1. / 3.));
    double rho_o_r = sqrt(2 * mu_r * Xo_r) * (Ri_r / hc);

    std::stringstream name;
    name << "Ed_BW_Er_" << E_BW << "_Gam_" << W_BW;

    auto Ed_BW = new TF1(
        name.str().c_str(), "([0]*sqrt([1]*x))/(pow((x-[2]),2) + pow(([3]*sqrt([1]*x)),2))", 0., fMax_energy_res);
    Ed_BW->SetParameter(0, Wo_r * (Ri_r / hc) / rho_o_r * Wo_r / 4);
    Ed_BW->SetParameter(1, 2. * mu_r);
    Ed_BW->SetParameter(2, Xo_r);
    Ed_BW->SetParameter(3, Wo_r * (Ri_r / hc) / rho_o_r / 2.);

    // For simulation with Breit Wigner in Sequential Decay
    double Xo_s = Er_BW; // Centroid BW
    double Wo_s = Wr_BW; // GAMMA BW
    // Breit Wigner parameters
    double mu_s = 931.5 * A / (A + 1);
    double Ri_s = 1.4 * (pow(A, 1. / 3.) + pow(1, 1. / 3.));
    double rho_o_s = sqrt(2 * mu_s * Xo_s) * (Ri_s / hc);

    std::stringstream name_s;
    name_s << "ESeq_BW_Er_" << Er_BW << "_Gam_" << Wr_BW;

    auto ESeq_BW = new TF1(
        name_s.str().c_str(), "([0]*sqrt([1]*x))/(pow((x-[2]),2) + pow(([3]*sqrt([1]*x)),2))", 0., fMax_energy_res);
    ESeq_BW->SetParameter(0, Wo_s * (Ri_s / hc) / rho_o_s * Wo_s / 4);
    ESeq_BW->SetParameter(1, 2. * mu_s);
    ESeq_BW->SetParameter(2, Xo_s);
    ESeq_BW->SetParameter(3, Wo_s * (Ri_s / hc) / rho_o_s / 2.);

    // For dineutron decay
    double Xo = E_DIN; // 0.1; //Centroid BW
    double Wo = W_DIN; // 0.1; //GAMMA BW

    // Breit Wigner parameters
    double mu = 931.5 * A / (A + 1);
    double Ri = 1.4 * (pow(A, 1. / 3.) + pow(1, 1. / 3.));
    double rho_o = sqrt(2 * mu * Xo) * (Ri / hc);

    auto BW_Dineutron =
        new TF1("BW_Dineutron", "([0]*sqrt([1]*x))/(pow((x-[2]),2) + pow(([3]*sqrt([1]*x)),2))", 0., 15.);
    BW_Dineutron->SetParameter(0, Wo * (Ri / hc) / rho_o * Wo / 4);
    BW_Dineutron->SetParameter(1, 2. * mu);
    BW_Dineutron->SetParameter(2, Xo);
    BW_Dineutron->SetParameter(3, Wo * (Ri / hc) / rho_o / 2.);

    auto tree = std::make_unique<TTree>("evt", "NeutronDecayTree");
    tree->Branch("Ed", &Ed);
    tree->Branch("Beta", &b_beam);

    if (N == 1)
    {
        tree->Branch("Pn_X", &Pn_X);
        tree->Branch("Pn_Y", &Pn_Y);
        tree->Branch("Pn_Z", &Pn_Z);
        tree->Branch("Pn_E", &Pn_E);
        tree->Branch("Pf_X", &Pf_X);
        tree->Branch("Pf_Y", &Pf_Y);
        tree->Branch("Pf_Z", &Pf_Z);
        tree->Branch("Pf_E", &Pf_E);
        tree->Branch("Estar_1n", &Estar_1n);
    }

    if (N == 2)
    {
        tree->Branch("m_nn", &m_nn);
        tree->Branch("m_fn1", &m_fn1); // with 1st neutron
        tree->Branch("m_fn2", &m_fn2); // with 2nd neutron
        tree->Branch("E_star2n", &E_star2n);

        tree->Branch("theta_n_n", &theta_n_n);
        tree->Branch("theta_fn2_n", &theta_fn2_n);
        tree->Branch("theta_fn1_n", &theta_fn1_n);
        tree->Branch("theta_nn1_f", &theta_nn1_f);
        tree->Branch("theta_nn2_f", &theta_nn2_f);
        tree->Branch("theta_kY_1", &theta_kY_1);
        tree->Branch("theta_kT_1", &theta_kT_1);
        tree->Branch("Ex_Y_1", &Ex_Y_1);
        tree->Branch("Ex_T_1", &Ex_T_1);
        tree->Branch("theta_kY_2", &theta_kY_2);
        tree->Branch("theta_kT_2", &theta_kT_2);
        tree->Branch("Ex_Y_2", &Ex_Y_2);
        tree->Branch("Ex_T_2", &Ex_T_2);
        tree->Branch("E_star2n_n1", &E_star2n_n1);
        tree->Branch("E_star2n_n2", &E_star2n_n2);
    }

    // Variables for C_nn test
    double Ran, Cnn_val;
    int Er_num = 1;
    int Gam_num = 1;

    int i_counter = 0;

    double Delta = 2.;

    // Nbr of simulated events loop
    for (Long64_t i = 0; i < Nb; i++)
    {
        showProgressBar(i, Nb);

        m_nn.clear();
        theta_n_n.clear();
        m_fn1.clear();
        m_fn2.clear();
        E_star2n.clear();

        E_star2n_n1.clear();
        E_star2n_n2.clear();
        theta_fn2_n.clear();
        theta_fn1_n.clear();
        theta_nn1_f.clear();
        theta_nn2_f.clear();
        theta_kY_1.clear();
        theta_kT_1.clear();
        Ex_Y_1.clear();
        Ex_T_1.clear();
        theta_kY_2.clear();
        theta_kT_2.clear();
        Ex_Y_2.clear();
        Ex_T_2.clear();

        r0 = r0_init;

        // Ed = 0.5;
        // Ed = Rand->Uniform(0.,12.);
        // Ed = Rand->Gaus(0.5,0.2);
        if (W_BW == 0.)
        {
            Ed = E_BW;
        }
        else if (W_BW == -1.)
        {
            Ed = Rand->Uniform(0., E_BW);
        }
        else if (W_BW == -2.)
        {
            Ed = Delta;
        }
        else
        {
            Ed = Ed_BW->GetRandom(0., fMax_energy_res);
        }

        if (i % 5000 == 0)
        {
            if (Delta < 1.5)
            {
                Delta = Delta + 0.5;
            }
            else
            {
                Delta = Delta + 1.;
            }
        }

        // Beta beam
        b_beam = Rand->Gaus(Beta, Beta_sig); //
        gamma_beam = 1. / sqrt(1. - b_beam * b_beam);

        /********************1 NEUTRON DECAY SIMULATION*********************/
        if (N == 1)
        {
            PHASE_2(mn, mf, Ed, Pn, Pf, Rand_angle.get());

            // Back in the Lab frame
            ZBOOST(Pn, b_beam);
            ZBOOST(Pf, b_beam);

            Pn_X = Pn[0];
            Pn_Y = Pn[1];
            Pn_Z = Pn[2];
            Pn_E = Pn[3];

            Pf_X = Pf[0];
            Pf_Y = Pf[1];
            Pf_Z = Pf[2];
            Pf_E = Pf[3];

            Estar_1n = Estar1n(mn, mf, Pn, Pf);

            // Generate Target Position
            double x = 0.;
            double y = 0.;
            double z = 0.;

            // double x = Rand->Gaus(0., 1.);       // cm
            // double y = Rand->Gaus(0., 1.);       // cm
            // double z = Rand->Uniform(-2.5, 2.5); // cm

            if (!isnan(Pf[2]))
            {
                outfile << i_counter << "  " << 2 << "\n";
                outfile << -1 << "  " << Z << "  " << A << "  " << Pf[0] / 1000. << "  " << Pf[1] / 1000. << "  "
                        << Pf[2] / 1000. << "  " << x << "  " << y << "  " << z << "\n";
                outfile << 2112 << "  " << 0 << "  " << 1 << "  " << Pn[0] / 1000. << "  " << Pn[1] / 1000. << "  "
                        << Pn[2] / 1000. << "  " << x << "  " << y << "  " << z << "\n";

                i_counter++;
            }
            tree->Fill();
        }

        /********************2 NEUTRONS DECAY SIMULATION*********************/
        if (N == 2)
        {

            t0 = h_barC_MeV_fm / Wr_BW;

            Cnn_max = Cnn_0(0., r0); // Normalization

            if (FSI_nn == 1)
            {
                do
                {
                    if (SEQ_Decay == 1)
                    { // Sequential decay mode
                        E23 = PHASE2x2(mn, mn, mf, Ed, Pnn1, Pnn2, Pnn3, Rand_angle.get(), A, ESeq_BW);

                        h1_E_res->Fill(E23);
                        q = sqrt(((Pnn1[0] - Pnn2[0]) * (Pnn1[0] - Pnn2[0]) +
                                  (Pnn1[1] - Pnn2[1]) * (Pnn1[1] - Pnn2[1]) +
                                  (Pnn1[2] - Pnn2[2]) * (Pnn1[2] - Pnn2[2])) -
                                 (Pnn1[3] - Pnn2[3]) * (Pnn1[3] - Pnn2[3])) /
                            2.;

                        Ran = Rand_0->Uniform(0., 1.);

                        C_nn_seq = Cnn_sequential(Pnn1, Pnn2, r0, t0);

                        // Cnn_val = Cnn_0(q,r0)/Cnn_max;

                        Cnn_val = C_nn_seq / Cnn_max;

                        h1_C_nn_seq->Fill(C_nn_seq);
                    }
                    else if (SEQ_Decay == 2)
                    { // Dineutron Decay

                        E23 = PHASE2x2(mf, mn, mn, Ed, Pnn3, Pnn2, Pnn1, Rand_angle.get(), A, BW_Dineutron);

                        q = sqrt(((Pnn1[0] - Pnn2[0]) * (Pnn1[0] - Pnn2[0]) +
                                  (Pnn1[1] - Pnn2[1]) * (Pnn1[1] - Pnn2[1]) +
                                  (Pnn1[2] - Pnn2[2]) * (Pnn1[2] - Pnn2[2])) -
                                 (Pnn1[3] - Pnn2[3]) * (Pnn1[3] - Pnn2[3])) /
                            2.;

                        Ran = 0.5;
                        Cnn_val = 1.;
                    }
                    else
                    { // Direct decay
                        PHASE_3(mn, mn, mf, Ed, Pnn1, Pnn2, Pnn3,
                                Rand_angle.get()); // direct decay

                        q = sqrt(((Pnn1[0] - Pnn2[0]) * (Pnn1[0] - Pnn2[0]) +
                                  (Pnn1[1] - Pnn2[1]) * (Pnn1[1] - Pnn2[1]) +
                                  (Pnn1[2] - Pnn2[2]) * (Pnn1[2] - Pnn2[2])) -
                                 (Pnn1[3] - Pnn2[3]) * (Pnn1[3] - Pnn2[3])) /
                            2.;

                        Ran = Rand_0->Uniform(0., 1.);

                        // Cnn_max = Cnn_0(0,r0);
                        Cnn_val = Cnn_0(q, r0) / Cnn_max;

                        h1_C_nn->Fill(Cnn_0(q, r0));
                    }
                } while (Ran > Cnn_val);
            }
            else
            {
                PHASE_3(mn, mn, mf, Ed, Pnn1, Pnn2, Pnn3, Rand_angle.get());
                q = sqrt(((Pnn1[0] - Pnn2[0]) * (Pnn1[0] - Pnn2[0]) + (Pnn1[1] - Pnn2[1]) * (Pnn1[1] - Pnn2[1]) +
                          (Pnn1[2] - Pnn2[2]) * (Pnn1[2] - Pnn2[2])) -
                         (Pnn1[3] - Pnn2[3]) * (Pnn1[3] - Pnn2[3])) /
                    2.;
            }

            // Back in the Lab frame
            ZBOOST(Pnn1, b_beam);
            ZBOOST(Pnn2, b_beam);
            ZBOOST(Pnn3, b_beam);

            // Making input file for R3BRoot simulation
            // Generate Target Position
            double x = 0.;
            double y = 0.;
            double z = 0.;

            // double x = Rand->Gaus(0., 1.);       // cm
            // double y = Rand->Gaus(0., 1.);       // cm
            // double z = Rand->Uniform(-2.5, 2.5); // cm

            outfile << i_counter << "  " << 3 << "\n";
            outfile << -1 << "  " << Z << "  " << A << "  " << Pnn3[0] / 1000. << "  " << Pnn3[1] / 1000. << "  "
                    << Pnn3[2] / 1000. << "  " << x << "  " << y << "  " << z << "\n";
            outfile << 2112 << "  " << 0 << "  " << 1 << "  " << Pnn1[0] / 1000. << "  " << Pnn1[1] / 1000. << "  "
                    << Pnn1[2] / 1000. << "  " << x << "  " << y << "  " << z << "\n";
            outfile << 2112 << "  " << 0 << "  " << 1 << "  " << Pnn2[0] / 1000. << "  " << Pnn2[1] / 1000. << "  "
                    << Pnn2[2] / 1000. << "  " << x << "  " << y << "  " << z << "\n";

            i_counter++;

            beta_n[0] = 1. / sqrt(1. + (mn * mn /
                                        (Pnn1[3] * Pnn1[3] - mn * mn - Pnn1[0] * Pnn1[0] -
                                         Pnn1[1] * Pnn1[1]))); // Beta_z 1st neutron
            beta_n[1] = 1. / sqrt(1. + (mn * mn /
                                        (Pnn2[3] * Pnn2[3] - mn * mn - Pnn2[0] * Pnn2[0] -
                                         Pnn2[1] * Pnn2[1]))); // Beta_z 2nd neutron

            gamma_n[0] = 1. / sqrt(1 - beta_n[0] * beta_n[0]);
            gamma_n[1] = 1. / sqrt(1 - beta_n[1] * beta_n[1]);

            // angle_n_n in CoM calculation
            P_n1[0] = Pnn1[0];
            P_n1[1] = Pnn1[1];
            P_n1[2] = Pnn1[2];

            P_n1[3] = TMath::Sqrt(mn * mn + P_n1[0] * P_n1[0] + P_n1[1] * P_n1[1] + P_n1[2] * P_n1[2]);

            P_n2[0] = Pnn2[0];
            P_n2[1] = Pnn2[1];
            P_n2[2] = Pnn2[2];

            P_n2[3] = TMath::Sqrt(mn * mn + P_n2[0] * P_n2[0] + P_n2[1] * P_n2[1] + P_n2[2] * P_n2[2]);

            P_f[0] = Pnn3[0];
            P_f[1] = Pnn3[1];
            P_f[2] = Pnn3[2];

            P_f[3] = TMath::Sqrt(mf * mf + P_f[0] * P_f[0] + P_f[1] * P_f[1] + P_f[2] * P_f[2]);

            P_tot[0] = P_n1[0] + P_n2[0] + P_f[0];
            P_tot[1] = P_n1[1] + P_n2[1] + P_f[1];
            P_tot[2] = P_n1[2] + P_n2[2] + P_f[2];

            P_tot[3] = TMath::Sqrt((mf + 2 * mn) * (mf + 2 * mn) + P_tot[0] * P_tot[0] + P_tot[1] * P_tot[1] +
                                   P_tot[2] * P_tot[2]);

            b_boost[0] = -P_tot[0] / P_tot[3];
            b_boost[1] = -P_tot[1] / P_tot[3];
            b_boost[2] = -P_tot[2] / P_tot[3];

            BOOST(P_n1, b_boost);
            BOOST(P_n2, b_boost);
            BOOST(P_f, b_boost);

            double theta_n_n_bfLAND =
                TMath::ACos(VoV(P_n1, P_n2) / (TMath::Sqrt(P_n1[0] * P_n1[0] + P_n1[1] * P_n1[1] + P_n1[2] * P_n1[2]) *
                                               TMath::Sqrt(P_n2[0] * P_n2[0] + P_n2[1] * P_n2[1] + P_n2[2] * P_n2[2])));

            P_n1n2[0] = -P_n1[0] + P_n2[0];
            P_n1n2[1] = -P_n1[1] + P_n2[1];
            P_n1n2[2] = -P_n1[2] + P_n2[2];

            P_n2n1[0] = -P_n1n2[0];
            P_n2n1[1] = -P_n1n2[1];
            P_n2n1[2] = -P_n1n2[2];

            P_fn1[0] = -P_f[0] + P_n1[0];
            P_fn1[1] = -P_f[1] + P_n1[1];
            P_fn1[2] = -P_f[2] + P_n1[2];

            P_fn2[0] = -P_f[0] + P_n2[0];
            P_fn2[1] = -P_f[1] + P_n2[1];
            P_fn2[2] = -P_f[2] + P_n2[2];

            double theta_fn1_n_bfLAND = TMath::ACos(
                VoV(P_fn1, P_n2) / (TMath::Sqrt(P_fn1[0] * P_fn1[0] + P_fn1[1] * P_fn1[1] + P_fn1[2] * P_fn1[2]) *
                                    TMath::Sqrt(P_n2[0] * P_n2[0] + P_n2[1] * P_n2[1] + P_n2[2] * P_n2[2])));
            double theta_fn2_n_bfLAND = TMath::ACos(
                VoV(P_n1, P_fn2) / (TMath::Sqrt(P_n1[0] * P_n1[0] + P_n1[1] * P_n1[1] + P_n1[2] * P_n1[2]) *
                                    TMath::Sqrt(P_fn2[0] * P_fn2[0] + P_fn2[1] * P_fn2[1] + P_fn2[2] * P_fn2[2])));
            double theta_nn1_f_bfLAND = TMath::ACos(
                VoV(P_n1n2, P_f) / (TMath::Sqrt(P_n1n2[0] * P_n1n2[0] + P_n1n2[1] * P_n1n2[1] + P_n1n2[2] * P_n1n2[2]) *
                                    TMath::Sqrt(P_f[0] * P_f[0] + P_f[1] * P_f[1] + P_f[2] * P_f[2])));
            double theta_nn2_f_bfLAND =
                TMath::ACos(VoV(P_f, P_n2n1) /
                            (TMath::Sqrt(P_f[0] * P_f[0] + P_f[1] * P_f[1] + P_f[2] * P_f[2]) *
                             TMath::Sqrt(P_n2n1[0] * P_n2n1[0] + P_n2n1[1] * P_n2n1[1] + P_n2n1[2] * P_n2n1[2])));

            kx_Y_1[0] = (mf * P_n1[0] - mn * P_f[0]) / (mn + mf);
            kx_Y_1[1] = (mf * P_n1[1] - mn * P_f[1]) / (mn + mf);
            kx_Y_1[2] = (mf * P_n1[2] - mn * P_f[2]) / (mn + mf);

            kx_T_1[0] = (mn * P_n1[0] - mn * P_n2[0]) / (mn + mn);
            kx_T_1[1] = (mn * P_n1[1] - mn * P_n2[1]) / (mn + mn);
            kx_T_1[2] = (mn * P_n1[2] - mn * P_n2[2]) / (mn + mn);

            ky_Y_1[0] = (mn * (P_f[0] + P_n1[0]) - (mn + mf) * P_n2[0]) / (2 * mn + mf);
            ky_Y_1[1] = (mn * (P_f[1] + P_n1[1]) - (mn + mf) * P_n2[1]) / (2 * mn + mf);
            ky_Y_1[2] = (mn * (P_f[2] + P_n1[2]) - (mn + mf) * P_n2[2]) / (2 * mn + mf);

            ky_T_1[0] = (mf * (P_n1[0] + P_n2[0]) - (mn + mn) * P_f[0]) / (2 * mn + mf);
            ky_T_1[1] = (mf * (P_n1[1] + P_n2[1]) - (mn + mn) * P_f[1]) / (2 * mn + mf);
            ky_T_1[2] = (mf * (P_n1[2] + P_n2[2]) - (mn + mn) * P_f[2]) / (2 * mn + mf);

            double theta_kY_1_bfLAND = TMath::ACos(VoV(kx_Y_1, ky_Y_1) / (VMOD(kx_Y_1) * VMOD(ky_Y_1)));
            double theta_kT_1_bfLAND = TMath::ACos(VoV(kx_T_1, ky_T_1) / (VMOD(kx_T_1) * VMOD(ky_T_1)));

            double Ex_T_1_bfLAND =
                (mn + mn) * (kx_T_1[0] * kx_T_1[0] + kx_T_1[1] * kx_T_1[1] + kx_T_1[2] * kx_T_1[2]) / (2 * mn * mn);
            double Ex_Y_1_bfLAND =
                (mn + mf) * (kx_Y_1[0] * kx_Y_1[0] + kx_Y_1[1] * kx_Y_1[1] + kx_Y_1[2] * kx_Y_1[2]) / (2 * mn * mf);

            kx_Y_2[0] = (mf * P_n2[0] - mn * P_f[0]) / (mn + mf);
            kx_Y_2[1] = (mf * P_n2[1] - mn * P_f[1]) / (mn + mf);
            kx_Y_2[2] = (mf * P_n2[2] - mn * P_f[2]) / (mn + mf);

            kx_T_2[0] = (mn * P_n2[0] - mn * P_n1[0]) / (mn + mn);
            kx_T_2[1] = (mn * P_n2[1] - mn * P_n1[1]) / (mn + mn);
            kx_T_2[2] = (mn * P_n2[2] - mn * P_n1[2]) / (mn + mn);

            ky_Y_2[0] = (mn * (P_f[0] + P_n2[0]) - (mn + mf) * P_n1[0]) / (2 * mn + mf);
            ky_Y_2[1] = (mn * (P_f[1] + P_n2[1]) - (mn + mf) * P_n1[1]) / (2 * mn + mf);
            ky_Y_2[2] = (mn * (P_f[2] + P_n2[2]) - (mn + mf) * P_n1[2]) / (2 * mn + mf);

            ky_T_2[0] = (mf * (P_n2[0] + P_n1[0]) - (mn + mn) * P_f[0]) / (2 * mn + mf);
            ky_T_2[1] = (mf * (P_n2[1] + P_n1[1]) - (mn + mn) * P_f[1]) / (2 * mn + mf);
            ky_T_2[2] = (mf * (P_n2[2] + P_n1[2]) - (mn + mn) * P_f[2]) / (2 * mn + mf);

            double theta_kY_2_bfLAND = TMath::ACos(VoV(kx_Y_2, ky_Y_2) / (VMOD(kx_Y_2) * VMOD(ky_Y_2)));
            double theta_kT_2_bfLAND = TMath::ACos(VoV(kx_T_2, ky_T_2) / (VMOD(kx_T_2) * VMOD(ky_T_2)));

            double Ex_T_2_bfLAND =
                (mn + mn) * (kx_T_2[0] * kx_T_2[0] + kx_T_2[1] * kx_T_2[1] + kx_T_2[2] * kx_T_2[2]) / (2 * mn * mn);
            double Ex_Y_2_bfLAND =
                (mn + mf) * (kx_Y_2[0] * kx_Y_2[0] + kx_Y_2[1] * kx_Y_2[1] + kx_Y_2[2] * kx_Y_2[2]) / (2 * mn * mf);

            E_star2n.push_back(Estar2n(mn, mf, Pnn1, Pnn2, Pnn3));
            theta_n_n.push_back(theta_n_n_bfLAND);
            theta_fn1_n.push_back(theta_fn1_n_bfLAND);
            theta_fn2_n.push_back(theta_fn2_n_bfLAND);
            theta_nn1_f.push_back(theta_nn1_f_bfLAND);
            theta_nn2_f.push_back(theta_nn2_f_bfLAND);
            theta_kY_1.push_back(theta_kY_1_bfLAND);
            theta_kT_1.push_back(theta_kT_1_bfLAND);
            Ex_T_1.push_back(Ex_T_1_bfLAND / E_star2n.at(0));
            Ex_Y_1.push_back(Ex_Y_1_bfLAND / E_star2n.at(0));
            theta_kY_2.push_back(theta_kY_2_bfLAND);
            theta_kT_2.push_back(theta_kT_2_bfLAND);
            Ex_T_2.push_back(Ex_T_2_bfLAND / E_star2n.at(0));
            Ex_Y_2.push_back(Ex_Y_2_bfLAND / E_star2n.at(0));
            E_star2n_n1.push_back(Estar1n(mn, mf, Pnn1, Pnn3));
            E_star2n_n2.push_back(Estar1n(mn, mf, Pnn2, Pnn3));
            m_nn.push_back(DALITZ(mn, mn, E_star2n.at(0), Pnn2, Pnn1));
            m_fn1.push_back(DALITZ(mf, mn, E_star2n.at(0), Pnn1, Pnn3));
            m_fn2.push_back(DALITZ(mf, mn, E_star2n.at(0), Pnn3, Pnn2));

            if (m_nn.size() != m_fn1.size() || m_nn.size() != m_fn2.size() || m_nn.size() != E_star2n.size())
            {
                std::cout << "Problem with the size of the vectors" << std::endl;
            }
            tree->Fill();
        }
    }

    // Store data in the ROOT file and close it
    TFile hfile(Output_Name + (TString) ".root", "RECREATE");
    tree->Write();
    hfile.Close();

    // Close ASCII file
    outfile.close();

    std::cout << std::endl;
    std::cout << "Calculation finished successfully" << std::endl << std::endl;
    gApplication->Terminate();
}
