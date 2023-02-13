#ifndef SampleHandler_h
#define SampleHandler_h

#include <TGraph.h>
#include <TGraph2DErrors.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

// clang-format off
#include "commons.h"
// clang-format on
#include "CepGen/Physics/Utils.h"

class SampleHandler {
public:
  static SampleHandler fromCLAS(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() < 5)
        continue;
      // q2, x, f2, stat, syst, total;
      hnd.addValue(dt.at(1), dt.at(0), dt.at(2), std::hypot(dt.at(3), dt.at(4)));
    }
    return hnd;
  }
  static SampleHandler fromHermes(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() != 8)
        continue;
      // x, q2, f2, stat, pid, model, mis, rad;
      const double err = sqrt(dt.at(3) * dt.at(3) + dt.at(4) * dt.at(4) + dt.at(5) * dt.at(5) + dt.at(6) * dt.at(6) +
                              dt.at(7) * dt.at(7));
      hnd.addValue(dt.at(0), dt.at(1), dt.at(2), err * dt.at(2) / 100.);
    }
    return hnd;
  }
  static SampleHandler fromZEUS(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() != 6)
        continue;
      // x, q2, f2, stat, sysp, sysm;
      hnd.addValue(dt.at(0), dt.at(1), dt.at(2), std::hypot(dt.at(3), 0.5 * (dt.at(4) + dt.at(5))));
    }
    return hnd;
  }
  static SampleHandler fromH1lowQ2(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() != 6)
        continue;
      // q2, x, f2, stat, syst, total;
      //const double err = std::hypot(stat, syst);
      hnd.addValue(dt.at(1), dt.at(0), dt.at(2), dt.at(2) * dt.at(5) / 100.);
    }
    return hnd;
  }
  static SampleHandler fromH11997(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() < 4)
        continue;
      // q2, x, f2, total (percent);
      hnd.addValue(dt.at(1), dt.at(0), dt.at(2), dt.at(2) * dt.at(3) / 100.);
    }
    return hnd;
  }
  static SampleHandler fromBCDMS(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() != 4)
        continue;
      // x, q2, f2, total;
      hnd.addValue(dt.at(0), dt.at(1), dt.at(2), dt.at(3));
    }
    return hnd;
  }
  static SampleHandler fromNMC(const std::string& filename) {
    SampleHandler hnd(filename);
    for (const auto& dt : hnd.raw_values_) {
      if (dt.size() != 5)
        continue;
      // x, q2, f2, stat, syst
      hnd.addValue(dt.at(0), dt.at(1), dt.at(2), std::hypot(dt.at(3), dt.at(4)));
    }
    return hnd;
  }

  SampleHandler& operator+(const SampleHandler& oth) {
    for (const auto& val : oth.values_)
      addValue(std::get<0>(val.first), std::get<1>(val.first), val.second.value, val.second.uncertainty);
    return *this;
  }

  void addValue(double xbj, double q2, double value, double unc = 0.) {
    values_[std::make_tuple(xbj, q2)] = value_t{value, unc};
  }
  TGraph2DErrors graph() const {
    TGraph2DErrors gr;
    size_t i = 0;
    for (const auto& val : values_) {
      gr.SetPoint(i, std::get<1>(val.first), std::get<0>(val.first), val.second.value);
      gr.SetPointError(i, 0., 0., val.second.uncertainty);
      ++i;
    }
    return gr;
  }
  TGraphErrors xBjGraph(double q2) const {
    TGraphErrors out;
    size_t i = 0;
    for (const auto& val : values_) {
      if (fabs(std::get<1>(val.first) - q2) / q2 > rel_diff_)
        continue;
      out.SetPoint(i, std::get<0>(val.first), val.second.value);
      out.SetPointError(i, 0., val.second.uncertainty);
      ++i;
    }
    return out;
  }
  TGraphErrors wGraph(double q2) const {
    TGraphErrors out;
    size_t i = 0;
    for (const auto& val : values_) {
      const auto &aQ2 = std::get<1>(val.first), &axBj = std::get<0>(val.first);
      if (fabs(aQ2 - q2) / q2 > rel_diff_)
        continue;
      const auto w = std::sqrt(cepgen::utils::mX2(axBj, aQ2, std::pow(0.938272, 2)));
      out.SetPoint(i, w, val.second.value);
      out.SetPointError(i, 0., val.second.uncertainty);
      ++i;
    }
    return out;
  }
  std::vector<double> q2Values() const {
    std::set<double> vals;
    for (const auto& val : values_)
      vals.insert(std::get<1>(val.first));
    return std::vector<double>(vals.begin(), vals.end());
  }

  TGraph wQ2Span() const {
    TGraph gr;
    for (const auto& val : values_) {
      const auto &q2 = std::get<1>(val.first), &xbj = std::get<0>(val.first);
      gr.AddPoint(std::sqrt(cepgen::utils::mX2(xbj, q2, std::pow(0.938272, 2))), q2);
    }
    return gr;
  }
  TGraph xbjQ2Span() const {
    TGraph gr;
    for (const auto& val : values_)
      gr.AddPoint(std::get<0>(val.first), std::get<1>(val.first));
    return gr;
  }

  TGraph2D wQ2Values() const {
    TGraph2D gr;
    for (const auto& val : values_) {
      const auto &q2 = std::get<1>(val.first), &xbj = std::get<0>(val.first);
      gr.AddPoint(std::sqrt(cepgen::utils::mX2(xbj, q2, std::pow(0.938272, 2))), q2, val.second.value);
    }
    return gr;
  }
  TGraph2D xbjQ2Values() const {
    TGraph2D gr;
    for (const auto& val : values_)
      gr.AddPoint(std::get<0>(val.first), std::get<1>(val.first), val.second.value);
    return gr;
  }

private:
  explicit SampleHandler(const std::string& filename) : filename_(filename) {
    std::ifstream cl(filename);
    if (!cl.is_open())
      throw std::runtime_error("Failed to open file '" + filename + "' for reading.");
    std::string line, word;
    while (std::getline(cl, line)) {
      if (line[0] == '#')
        continue;
      auto& raw_line = raw_values_.emplace_back();
      std::stringstream ss(line);
      while (std::getline(ss, word, ','))
        raw_line.emplace_back(std::stod(word));
    }
  }
  struct value_t {
    double value, uncertainty;
  };
  const std::string filename_;
  std::vector<std::vector<double> > raw_values_;
  std::map<std::tuple<double, double>, value_t> values_;
  const double rel_diff_ = 0.01;
};

#endif
