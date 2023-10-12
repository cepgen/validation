// clang-format off
#include "CepGenEnvironment.h"
R__LOAD_LIBRARY(libCepGenLHAPDF.so)
R__LOAD_LIBRARY(libCepGenPythia8.so)
// clang-format on

#include "CepGen/CollinearFluxes/CollinearFlux.h"
#include "CepGen/Modules/PartonFluxFactory.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace cepgen;

void coll_fluxes_2d(vector<string> types = {},
                    bool log_q2 = true,
                    bool log_xbj = false,
                    const string& basename = "test") {
  pdg::MCDFileParser::parse("/home/laurent/work/dev/cepgen/External/mass_width_2023.txt");

  if (types.empty())
    types = CollinearFluxFactory::get().modules();

  //const Limits xbj_lim{1.e-5, 1.};
  //const Limits q2_lim{1.e-1, 100.};
  //const Limits xbj_lim{5.e-6, 0.9};
  //const Limits q2_lim{1., 200.};
  const Limits xbj_lim{5.e-4, 0.999}, q2_lim{0.1, 500.};  // hybrid SFs
  //const Limits xbj_lim{5.e-5, 0.9}, q2_lim{1., 500.};  // perturbative SFs
  //const Limits xbj_lim{5.e-4, 0.999}, q2_lim{0.1, 10.};
  //const Limits xbj_lim{5.e-2, 0.9}, q2_lim{5.e-2, 1.};  // resonances SFs
  //const Limits xbj_lim{5.e-3, 0.99}, q2_lim{5.e-2, 1.};
  const unsigned int num_points_q2 = 200, num_points_xbj = 200;

  for (const auto& type : types) {
    TGraph2D gr_flux;
    unsigned int idx = 0;
    const auto& flux = CollinearFluxFactory::get().build(type);
    for (const auto& q2 : q2_lim.generate(num_points_q2, log_q2)) {
      for (const auto& xbj : xbj_lim.generate(num_points_xbj, log_xbj)) {
        //cout << q2 << "|" << xbj << "|" << flux->fluxQ2(xbj, q2) << "\n";
        gr_flux.SetPoint(idx, log_xbj ? log10(xbj) : xbj, log_q2 ? log10(q2) : q2, flux->fluxQ2(xbj, q2));
        idx++;
      }
    }
    gStyle->SetPalette(kColorPrintableOnGrey);
    gStyle->SetPalette(kLightTemperature);
    TColor::InvertPalette();
    //gStyle->SetNdivisions(99, "Z");
    ostringstream os;
    os << CollinearFluxFactory::get().describe(type);

    ROOTCanvas c(basename + "_" + type, os.str());
    //c.SetLogz();
    gr_flux.Draw("surf3");
    string xbj_axis = "x_{Bj}", q2_axis = "Q^{2} (GeV^{2})";
    if (log_q2)
      q2_axis = "log_{10}(Q^{2}/GeV^{2})";
    if (log_xbj)
      xbj_axis = "log_{10}(" + xbj_axis + ")";
    gr_flux.GetHistogram()->SetTitle(Form(";%s;%s;f(x,Q^{2})", xbj_axis.data(), q2_axis.data()));
    //gr_flux.GetHistogram()->SetContour(90);
    c.Prettify(gr_flux.GetHistogram());
    gr_flux.GetXaxis()->SetTitleOffset(1.4);
    gr_flux.GetYaxis()->SetTitleOffset(1.8);
    gr_flux.GetZaxis()->SetTitleOffset(1.2);
    if (gr_flux.GetHistogram()->GetMaximum() > 1.)
      gr_flux.SetMaximum(1.);
    c.Save("png");
  }
}
