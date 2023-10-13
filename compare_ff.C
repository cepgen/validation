// clang-format off
#include "CepGenEnvironment.h"
// clang-format on

#include "CepGen/Core/Exception.h"
#include "CepGen/FormFactors/Parameterisation.h"
#include "CepGen/Generator.h"
#include "CepGen/Modules/FormFactorsFactory.h"
#include "CepGen/Physics/Beam.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"

void compare_ff(bool logx = false, bool logy = false, bool right = false, bool show_legend = true) {
  const cepgen::Limits q2_lims{1., 1.e5};
  const unsigned short num_points = 1000;

  pdg::MCDFileParser::parse(std::string(CEPGEN_PATH) + "/External/mass_width_2023.txt");

  cepgen::ROOTCanvas c("ff_comparison", "CepGen form factors");
  c.SetLegendX1(0.4);
  c.SetLegendY1(0.77);

  TMultiGraph mg;

  vector<int> colours = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1};
  vector<TGraph*> v_g_fe, v_g_fm;
  vector<pair<string, string> > v_p_ffnames = {{"StandardDipole", "Standard dipole"},
                                               {"Mergell", "Mergell et al."},
                                               {"Brash", "Brash et al."},
                                               {"Arrington", "Arrington et al."}};
  size_t j = 0;
  CG_INFO("") << v_p_ffnames;
  for (const auto& p_ffnames : v_p_ffnames) {
    const auto& ffmod = cepgen::FormFactorsFactory::get().build(p_ffnames.first);
    auto gr_fe = new TGraph();
    gr_fe->SetLineColor(colours[j]);
    gr_fe->SetLineWidth(2);
    gr_fe->SetLineStyle(1);
    mg.Add(gr_fe);
    v_g_fe.emplace_back(gr_fe);
    auto gr_fm = new TGraph();
    gr_fm->SetLineColor(colours[j]);
    gr_fm->SetLineWidth(2);
    gr_fm->SetLineStyle(2);
    mg.Add(gr_fm);
    v_g_fm.emplace_back(gr_fm);
    if (show_legend)
      c.AddLegendEntry(gr_fe, p_ffnames.second.c_str(), "l");
    for (const auto& q2 : q2_lims.generate(num_points, logx)) {
      auto fn = (*ffmod)(q2);
      gr_fe->SetPoint(gr_fe->GetN(), q2, fn.FE);
      gr_fm->SetPoint(gr_fm->GetN(), q2, fn.FM);
    }
    ++j;
  }

  //--- general plotting

  mg.Draw("al");
  mg.GetHistogram()->SetTitle("Q^{2}\\F_{E,M}^{p}");
  c.Prettify(mg.GetHistogram());
  mg.GetHistogram()->SetTitle("");

  mg.GetHistogram()->GetXaxis()->SetLimits(q2_lims.min(), q2_lims.max());
  if (logx)
    c.SetLogx();
  if (logy) {
    mg.GetHistogram()->SetMinimum(1.e-3);
    mg.GetHistogram()->SetMaximum(TMath::Min(5., mg.GetHistogram()->GetMaximum()));
    c.SetLogy();
  }
  //else mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, 1.5 );
  else {
    //mg.GetHistogram()->SetMinimum( 1.e-4 );
    mg.GetHistogram()->GetYaxis()->SetRangeUser(1.e-4, TMath::Min(1.5, mg.GetHistogram()->GetMaximum()));
    //mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, TMath::Min( 1.9, mg.GetHistogram()->GetMaximum() ) );
  }

  c.Save("pdf");
}
