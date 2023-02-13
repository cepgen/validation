// clang-format off
#include "commons.h"
// clang-format on

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TSystem.h>

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "SampleHandler.h"

void compare_sf_vs_w(double q2 = 1.225, bool logx = false, bool logy = false, bool plot_fl = false) {
  pdg::MCDFileParser::parse(cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mass_width_2021.mcd");

  const auto mp = cepgen::PDG::get().mass(2212);
  const cepgen::Limits w_lims{1., 1000.};
  const unsigned short num_points = 5000;
  const double x_pos_leg_data = 0.16, y_pos_leg_data = 0.825;

  const string mstw_grid_path = cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mstw_sf_scan_nnlo.dat";

  vector<StrFunParams> sfs;
  sfs.emplace_back(
      301,
      cepgen::ParametersList().set<cepgen::ParametersList>(
          "perturbativeSF",
          cepgen::ParametersList().setName<int>(205 /* MSTWgrid */).set<string>("gridPath", mstw_grid_path))
      //.set<vector<double> >("W2limits", {3., 16.})
  );
  sfs.emplace_back(11);
  sfs.emplace_back(101);
  sfs.emplace_back(103);
  sfs.emplace_back(102);
  sfs.emplace_back(202);
  sfs.emplace_back(203);
  sfs.emplace_back(204);
  //sfs.emplace_back(304);

  for (const auto& w : w_lims.generate(num_points, logx)) {
    const auto xbj = cepgen::utils::xBj(q2, mp * mp, w * w);
    for (auto& sf : sfs) {
      if (plot_fl)
        sf.graph().AddPoint(w, sf.sf()->FL(xbj, q2));
      else
        sf.graph().AddPoint(w, sf.sf()->F2(xbj, q2));
    }
  }

  cepgen::ROOTCanvas c("sf_comparison");
  c.SetTopLabel(Form("Q^{2} ~ %g GeV^{2}", q2));
  c.SetLegendX1(0.4);
  c.SetLegendY1(0.77);

  TMultiGraph mg;

  //--- model curves

  for (auto& sf : sfs) {
    if (sf.name().find("LUXlike") != std::string::npos) {
      auto gr_luxlike_line = (TGraph*)sf.graph().Clone();
      gr_luxlike_line->SetLineColorAlpha(kRed + 1, 0.5);
      gr_luxlike_line->SetLineWidth(8);
      mg.Add(gr_luxlike_line);
    }
    sf.graph().SetFillStyle(0);
    sf.graph().SetFillColor(0);
    mg.Add(&sf.graph());
    c.AddLegendEntry(&sf.graph(), cepgen::utils::split(sf.name(), ' ')[0].data(), "l");
  }

  if (sfs.size() > 9) {
    auto leg = c.GetLegend();
    leg->SetNColumns(3);
    leg->SetX1(0.18);
    leg->SetX2(0.88);
    leg->SetY1(0.77);
    leg->SetY2(0.92);
  }

  struct DataSample {
    explicit DataSample(const SampleHandler& hnd, const std::string& aname, int astyle)
        : handler(hnd), name(aname), style(astyle) {}
    SampleHandler handler;
    std::string name;
    int style{0};
  };
  std::vector<DataSample> samples;
  samples.emplace_back(SampleHandler::fromCLAS("samples/clas_0p225-4p725gev2.csv"), "CLAS", 24);
  samples.emplace_back(SampleHandler::fromBCDMS("samples/bcdms.csv"), "BCDMS", 30);
  samples.emplace_back(
      SampleHandler::fromCLAS("samples/nmc_0p75-65gev2.csv") + SampleHandler::fromNMC("samples/nmc_0p8-62gev2.csv"),
      "NMC",
      25);
  samples.emplace_back(SampleHandler::fromCLAS("samples/e665.csv"), "E665", 27);
  samples.emplace_back(SampleHandler::fromZEUS("samples/zeus.csv") +
                           SampleHandler::fromCLAS("samples/zeus_0p11-0p65gev2.csv") +
                           SampleHandler::fromCLAS("samples/zeus_1p5-15gev2.csv") +
                           SampleHandler::fromZEUS("samples/zeus_0p6-17gev2.csv") +
                           SampleHandler::fromCLAS("samples/zeus_8p5-5000gev2.csv"),
                       "ZEUS",
                       26);
  samples.emplace_back(SampleHandler::fromCLAS("samples/h1.csv") +
                           SampleHandler::fromCLAS("samples/h1_4p5-1600gev2.csv") +
                           //SampleHandler::fromCLAS("samples/h1_neutral_2004.csv") +
                           //SampleHandler::fromCLAS("samples/hera-h1-nc.csv") +
                           SampleHandler::fromH11997("samples/h1_1p5-150gev2.csv") +
                           SampleHandler::fromH1lowQ2("samples/h1_compton.csv") +
                           SampleHandler::fromH1lowQ2("samples/h1_0p35-3p5gev2.csv"),
                       "H1",
                       22);
  samples.emplace_back(SampleHandler::fromHermes("samples/hermes.csv"), "HERMES", 28);

  auto leg_data = new TLegend(x_pos_leg_data, y_pos_leg_data, x_pos_leg_data + 0.15, y_pos_leg_data + 0.1, "");
  if (!plot_fl)
    for (const auto& sample : samples) {
      auto* proj = new TGraphErrors(sample.handler.wGraph(q2));
      proj->SetMarkerStyle(sample.style);
      proj->SetMarkerSize(0.9);
      mg.Add(proj, "p");
      if (proj->GetN() > 0)
        leg_data->AddEntry(proj, sample.name.data(), "ep");
    }

  //--- general plotting

  mg.Draw("al");
  mg.GetHistogram()->SetTitle(Form("W (GeV)\\%s", (plot_fl ? "F_{L}^{p}" : "F_{2}^{p}")));
  c.Prettify(mg.GetHistogram());
  mg.GetHistogram()->SetTitle("");

  mg.GetHistogram()->GetXaxis()->SetLimits(w_lims.min(), w_lims.max());
  if (logx)
    c.SetLogx();
  if (logy) {
    mg.GetHistogram()->SetMinimum(1.e-3);
    c.SetLogy();
  } else {
    const double y_max = (plot_fl ? 5.0 : 1.5);
    mg.GetHistogram()->GetYaxis()->SetRangeUser(1.e-4, TMath::Min(y_max, mg.GetHistogram()->GetMaximum()));
  }
  if (leg_data->GetNRows() > 0) {
    leg_data->SetTextFont(132);
    leg_data->SetTextSize(0.035);
    leg_data->SetLineWidth(0);
    leg_data->SetFillStyle(0);
    leg_data->Draw();
  }
  c.SetGrid(1, 1);
  c.Save("pdf");
}
