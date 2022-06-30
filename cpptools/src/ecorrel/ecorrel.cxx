#include "ecorrel.hh"
#include <cmath>
#include <stdexcept>

namespace EnergyCorrelators
{
    CorrelatorsContainer::CorrelatorsContainer()
    : fr()
    , fw()
    , frxw()
    {
        ;
    }
    
    CorrelatorsContainer::~CorrelatorsContainer()
    {
        ;
    }

    void CorrelatorsContainer::clear()
    {
        fw.clear();
        fr.clear();
        frxw.clear();
    }

    void CorrelatorsContainer::addwr(const double &w, const double &r)
    {
        fw.push_back(w);
        fr.push_back(r);
    }

    std::vector<double> *CorrelatorsContainer::weights()
    {
        return &fw;
    }

    std::vector<double> *CorrelatorsContainer::rs()
    {
        return &fr;
    }

    const double *CorrelatorsContainer::wa()
    {
        return &fw[0];
    }

    const double *CorrelatorsContainer::ra()
    {
        return &fr[0];
    }

    std::vector<double> *CorrelatorsContainer::rxw()
    {
        frxw.clear();
        for (size_t i = 0; i < fr.size(); i++)
        {
            frxw.push_back(fr[i] * fw[i]);
        }
        return &frxw;
    }

    std::vector<fastjet::PseudoJet> constituents_as_vector(const fastjet::PseudoJet &jet)
    {
        std::vector<fastjet::PseudoJet> _v;
        for (auto &c : jet.constituents())
        {
            _v.push_back(c);
        }
        return _v;
    }

    CorrelatorBuilder::CorrelatorBuilder()
    : fec()
    , fncmax(4)
    {
        for (int i = 0; i < fncmax - 2 + 1; i++)
        {
            fec[i] = new CorrelatorsContainer();
        }
    }

    CorrelatorBuilder::CorrelatorBuilder(const std::vector<fastjet::PseudoJet> &parts, const double &scale, int nmax)
    : fec()
    , fncmax(nmax)
    {
        if (fncmax < 2)
        {
            throw std::overflow_error("asking for n-point correlator with n < 2?");
        }
        if (fncmax > 4)
        {
            throw std::overflow_error("max n for n-point correlator is currently 4");
        }
        for (int i = 0; i < fncmax - 2 + 1; i++)
        {
            fec[i] = new CorrelatorsContainer();
        }
        for (size_t i = 0; i < parts.size(); i++)
        {
            for (size_t j = i + 1; j < parts.size(); j++)
            {
                double _d12 = parts[i].delta_R(parts[j]);
                double _w2 = parts[i].E() * parts[j].E() / std::pow(scale, 2);
                fec[2 - 2]->addwr(_w2, _d12);
                if (fncmax < 3)
                    continue;
                for (size_t k = j + 1; k < parts.size(); k++)
                {
                    double _d13 = parts[i].delta_R(parts[k]);
                    double _d23 = parts[j].delta_R(parts[k]);
                    double _w3 = parts[i].E() * parts[j].E() * parts[k].E() / std::pow(scale, 3);
                    double _d3max = std::max(std::max(_d13, _d23), _d12);
                    fec[3 - 2]->addwr(_w3, _d3max);
                    if (fncmax < 4)
                        continue;
                    for (size_t l = k + 1; l < parts.size(); l++)
                    {
                        double _d14 = parts[i].delta_R(parts[l]);
                        double _d24 = parts[j].delta_R(parts[l]);
                        double _d34 = parts[k].delta_R(parts[l]);
                        double _w4 = parts[i].E() * parts[j].E() * parts[k].E() * parts[l].E() / std::pow(scale, 4);
                        double _d4max = std::max(std::max(std::max(_d13, _d23), std::max(_d14, _d24)), std::max(_d34, _d12));
                        fec[4 - 2]->addwr(_w4, _d4max);
                    }
                }
            }
        }
    }

    CorrelatorsContainer* CorrelatorBuilder::correlator(int n)
    {
        if (n > fncmax)
        {
            throw std::overflow_error("requesting n-point correlator with too large n");
        }
        if (n < 2)
        {
            throw std::overflow_error("requesting n-point correlator with n < 2?");
        }
        return fec[n - 2];
    }

    CorrelatorBuilder::~CorrelatorBuilder()
    {
        for (int i = 0; i < fncmax - 2 + 1; i++)
        {
            delete fec[i];
        }
    }

}
