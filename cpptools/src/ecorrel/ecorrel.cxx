#include "ecorrel.hh"
#include <cmath>

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

    CorrelatorsContainer EECrw(const std::vector<fastjet::PseudoJet> &parts, const double &scale)
    {   
        CorrelatorsContainer _c;
        for (size_t i = 0; i < parts.size(); i++)
        {
            for (size_t j = i+1; j < parts.size(); j++)
            {
                double _d12 = parts[i].delta_R(parts[j]);
                double _w = parts[i].E() * parts[j].E() / std::pow(scale, 2);
                double _r = _d12 * _w;
                _c.addwr(_w, _d12);
            }
        }
        return _c;
    }

    CorrelatorsContainer E3Crw(const std::vector<fastjet::PseudoJet> &parts, const double &scale)
    {
        CorrelatorsContainer _c;
        for (size_t i = 0; i < parts.size(); i++)
        {
            for (size_t j = i + 1; j < parts.size(); j++)
            {
                double _d12 = parts[i].delta_R(parts[j]);
                for (size_t k = j + 1; k < parts.size(); k++)
                {
                    double _d13 = parts[i].delta_R(parts[k]);
                    double _d23 = parts[j].delta_R(parts[k]);
                    double _w = parts[i].E() * parts[j].E() * parts[k].E() / std::pow(scale, 3);
                    double _d = std::max(std::max(_d13, _d23), _d12);
                    _c.addwr(_w, _d);
                }
            }
        }
        return _c;
    }

    CorrelatorsContainer E4Crw(const std::vector<fastjet::PseudoJet> &parts, const double &scale)
    {
        CorrelatorsContainer _c;
        for (size_t i = 0; i < parts.size(); i++)
        {
            for (size_t j = i + 1; j < parts.size(); j++)
            {
                double _d12 = parts[i].delta_R(parts[j]);
                for (size_t k = j + 1; k < parts.size(); k++)
                {
                    double _d13 = parts[i].delta_R(parts[k]);
                    double _d23 = parts[j].delta_R(parts[k]);
                    for (size_t l = k + 1; l < parts.size(); l++)
                    {
                        double _d14 = parts[i].delta_R(parts[l]);
                        double _d24 = parts[j].delta_R(parts[l]);
                        double _d34 = parts[k].delta_R(parts[l]);
                        double _w = parts[i].E() * parts[j].E() * parts[k].E() * parts[l].E() / std::pow(scale, 4);
                        double _d = std::max(std::max(std::max(_d13, _d23), std::max(_d14, _d24)), std::max(_d34, _d12));
                        _c.addwr(_w, _d);
                    }
                }
            }
        }
        return _c;
    }
}
