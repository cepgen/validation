#ifndef StrFunParams_h
#define StrFunParams_h

// clang-format off
#include "CepGenEnvironment.h"
// clang-format on

#include <TGraph.h>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>

#include <memory>

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

map<int, int> colours = {{301, kRed + 1},
                         {11, kCyan - 6},
                         {101, kGreen + 2},
                         {103, kBlue + 1},
                         {102, kMagenta + 1},
                         {202, kBlack},
                         {203, kBlack},
                         {204, kBlack},
                         {303, kOrange},
                         {401, kGray + 1},
                         {402, kTeal}};
map<int, int> styles = {
    {301, 1}, {11, 1}, {101, 1}, {103, 1}, {102, 1}, {202, 1}, {203, 3}, {204, 2}, {303, 1}, {401, 1}, {402, 1}};

class StrFunParams {
public:
  explicit StrFunParams(int type,
                        const cepgen::ParametersList& params = cepgen::ParametersList(),
                        const std::string& name = "",
                        int colour = -1,
                        int style = -1)
      : name_(name.empty() ? cepgen::StructureFunctionsFactory::get().describe(type) : name),
        sf_(cepgen::StructureFunctionsFactory::get().build(type, params)) {
    graph_.SetLineColor(colour >= 0 ? colour : colours.count(type) ? colours.at(type) : kBlack);
    graph_.SetLineStyle(style >= 0 ? style : styles.count(type) ? styles.at(type) : 1);
    graph_.SetLineWidth(3);
    graph_.SetFillStyle(0);
    graph_.SetFillColor(0);
  }

  const std::string& name() const { return name_; }
  cepgen::strfun::Parameterisation* sf() { return sf_.get(); }

  TGraph& graph() { return graph_; }

private:
  const std::string name_;
  const std::shared_ptr<cepgen::strfun::Parameterisation> sf_;
  TGraph graph_;
};

#endif
