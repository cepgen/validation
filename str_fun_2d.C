// clang-format off
#include "CepGenEnvironment.h"
R__LOAD_LIBRARY(libCepGenLHAPDF.so)
// clang-format on

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"

using namespace cepgen;

void str_fun_2d(int type, bool log_q2 = true, bool log_xbj = true, const char* basename = "test") {
  pdg::MCDFileParser::parse("/home/laurent/work/dev/cepgen/External/mass_width_2023.txt");

  TGraph2D gr_f2, gr_fl;
  //const Limits xbj_lim{1.e-5, 1.};
  //const Limits q2_lim{1.e-1, 100.};
  //const Limits xbj_lim{5.e-6, 0.9};
  //const Limits q2_lim{1., 200.};
  const Limits xbj_lim{5.e-4, 0.999}, q2_lim{0.1, 500.};  // hybrid SFs
  //const Limits xbj_lim{5.e-5, 0.9}, q2_lim{1., 500.}; // perturbative SFs
  //const Limits xbj_lim{5.e-4, 0.999}, q2_lim{0.1, 10.};
  //const Limits xbj_lim{5.e-2, 0.9}, q2_lim{5.e-2, 1.};// resonances SFs
  //const Limits xbj_lim{5.e-3, 0.99}, q2_lim{5.e-2, 1.};
  const unsigned int num_points_q2 = 200, num_points_xbj = 200;

  unsigned int idx = 0;
  const auto& str_fun = StructureFunctionsFactory::get().build(
      type, ParametersList().set<bool>("higherTwist", true).set<string>("pdfSet", "LUXlep-NNPDF31_nlo_as_0118_luxqed"));
  for (const auto& q2 : q2_lim.generate(num_points_q2, log_q2)) {
    for (const auto& xbj : xbj_lim.generate(num_points_xbj, log_xbj)) {
      const auto lxbj = log10(xbj), lq2 = log10(q2);
      cout << q2 << "|" << xbj << "|" << str_fun->operator()(xbj, q2) << "\n";
      gr_f2.SetPoint(idx, log_xbj ? lxbj : xbj, log_q2 ? lq2 : q2, str_fun->F2(xbj, q2));
      gr_fl.SetPoint(idx, log_xbj ? lxbj : xbj, log_q2 ? lq2 : q2, str_fun->FL(xbj, q2));
      idx++;
    }
  }

  //auto f_w2 = new TF1( "f_w2", "log10((10**x)/(1.-(10**x))*([0]-[1]**2))", lxbj_min, lxbj_max );
  /*auto f_w2 = new TF1("f_w2", "log10((10**y)/(1.-(10**y))*([0]-[1]**2))", lxbj_min, lxbj_max);
  f_w2->SetParameter(1, ParticleProperties::mass(PDG::Proton));
  f_w2->SetLineColor(kGreen + 1);
  f_w2->SetRange(lxbj_min, lxbj_max);*/
  /*auto f_q2 = new TF1("f_q2", "[0]", lxbj_min, lxbj_max);
  f_q2->SetParameter(0, log10(9.));
  //  f_q2->SetParameter( 0, 0.1 );
  f_q2->SetLineColor(kGray + 1);
  f_q2->SetRange(lxbj_min, lxbj_max);
  cout << f_q2->Eval(lxbj_min) << "|" << f_q2->Eval(lxbj_max) << endl;*/
  //TF1 f_q2("f_q2", "log10(9.)", lxbj_min, lxbj_max);
  gStyle->SetPalette(kColorPrintableOnGrey);
  gStyle->SetPalette(kLightTemperature);
  TColor::InvertPalette();
  //gStyle->SetNdivisions(99, "Z");
  std::ostringstream os;
  os << StructureFunctionsFactory::get().describe(type);
  string xbj_axis = "x_{Bj}", q2_axis = "Q^{2} (GeV^{2})";
  if (log_q2)
    q2_axis = "log_{10}(Q^{2}/GeV^{2})";
  if (log_xbj)
    xbj_axis = "log_{10}(" + xbj_axis + ")";
  {
    ROOTCanvas c(Form("%s_f2", basename), Form("%s F_{2}^{p} structure function", os.str().c_str()));
    //c.SetLogx();
    //c.SetLogy();
    gr_f2.Draw("surf3");
    gr_f2.GetHistogram()->SetTitle(Form(";%s;%s;F_{2}^{p}(x_{Bj},Q^{2})", xbj_axis.data(), q2_axis.data()));
    //gr_f2.GetHistogram()->SetContour(90);
    //gr_f2.SetMaximum(6.);
    /*f_w2->SetParameter(0, 3.);
    f_w2->DrawCopy("same");
    f_w2->SetParameter(0, 4.);
    f_w2->DrawCopy("same");
    f_q2.DrawCopy("same");*/
    c.Prettify(gr_f2.GetHistogram());
    gr_f2.GetXaxis()->SetTitleOffset(1.4);
    gr_f2.GetYaxis()->SetTitleOffset(1.8);
    gr_f2.GetZaxis()->SetTitleOffset(1.2);
    /*gr_f2.GetXaxis()->SetRangeUser(xbj_lim.min(), xbj_lim.max());
    gr_f2.GetYaxis()->SetRangeUser(q2_lim.min(), q2_lim.max());
    gPad->SetPhi(M_PI / 3);*/
    c.Save("png");
  }
  {
    ROOTCanvas c(Form("%s_fl", basename), Form("%s F_{L}^{p} structure function", os.str().c_str()));
    gr_fl.Draw("surf3");
    //    gr_fl.Draw( "lego2 0" );
    gr_fl.GetHistogram()->SetTitle(Form(";%s;%s;F_{L}^{p}(x_{Bj},Q^{2})", xbj_axis.data(), q2_axis.data()));
    //gr_fl.GetHistogram()->SetContour( 90 );
    //gr_fl.SetMaximum( 1.2 );
    //c.SetLogx();
    //c.SetLogy();
    //c.SetLogz();
    /*f_w2->SetParameter( 0, 3. );
    f_w2->DrawCopy( "same" );
    f_w2->SetParameter( 0, 4. );
    f_w2->DrawCopy( "same" );
    f_q2.DrawCopy( "same" );*/
    c.Prettify(gr_fl.GetHistogram());
    gr_fl.GetXaxis()->SetTitleOffset(1.4);
    gr_fl.GetYaxis()->SetTitleOffset(1.8);
    gr_fl.GetZaxis()->SetTitleOffset(1.2);
    c.Save("png");
  }
}
