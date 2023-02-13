// clang-format off
#include "CepGenEnvironment.h"
R__LOAD_LIBRARY(libCepGenLHAPDF.so)
// clang-format on

#include "Canvas.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/StructureFunctions/Parameterisation.h"

using namespace cepgen;

void f2_2d(int type, const char* basename = "test") {
  pdg::MCDFileParser::parse("/home/laurent/work/dev/cepgen/External/mass_width_2021.mcd");

  TGraph2D gr_f2, gr_fl;
  //const double lxbj_min = -5., lxbj_max = 0.;
  //const double lq2_min = -1., lq2_max = 2.;
  //const double xbj_min = 5.e-6, xbj_max = 0.9;
  //const double q2_min = 1., q2_max = 200.;
  const double xbj_min = 5.e-4, xbj_max = 0.999, q2_min = 0.1, q2_max = 500.;  // hybrid SFs
  //const double xbj_min = 5.e-5, xbj_max = 0.9, q2_min = 1., q2_max = 500.; // perturbative SFs
  //const double xbj_min = 5.e-4, xbj_max = 0.999, q2_min = 0.1, q2_max = 10.;
  //const double xbj_min = 5.e-2, xbj_max = 0.9, q2_min = 5.e-2, q2_max = 1.;// resonances SFs
  //const double xbj_min = 5.e-3, xbj_max = 0.99, q2_min = 5.e-2, q2_max = 1.;
  const double lxbj_min = log10(xbj_min), lxbj_max = log10(xbj_max);
  const double lq2_min = log10(q2_min), lq2_max = log10(q2_max);
  const unsigned int num_points_q2 = 200, num_points_xbj = 200;

  unsigned int idx = 0;
  const auto& str_fun = StructureFunctionsFactory::get().build(
      type, ParametersList().set<bool>("higherTwist", true).set<string>("pdfSet", "LUXlep-NNPDF31_nlo_as_0118_luxqed"));
  /*  auto str_fun = make_shared<strfun::Schaefer>();
  str_fun->params.perturbative_model.reset( new strfun::LHAPDF( "cteq6l1" ) );*/
  for (unsigned int i = 0; i < num_points_q2; ++i) {
    //const double q2 = pow( 10., lq2_min + i * ( lq2_max-lq2_min ) / ( num_points_q2-1 ) ), lq2 = log10( q2 );
    //const double q2 = q2_min + i * ( q2_max-q2_min ) / ( num_points_q2-1 );
    const double lq2 = lq2_min + i * (lq2_max - lq2_min) / (num_points_q2 - 1), q2 = pow(10., lq2);
    for (unsigned int j = 0; j < num_points_xbj; ++j) {
      //const double xbj = pow( 10., lxbj_min + j * ( lxbj_max-lxbj_min ) / ( num_points_xbj-1 ) );
      //const double xbj = xbj_min + j * ( xbj_max-xbj_min ) / ( num_points_xbj-1 );
      //const double lxbj = log10( xbj );
      const double lxbj = lxbj_min + j * (lxbj_max - lxbj_min) / (num_points_xbj - 1), xbj = pow(10., lxbj);
      auto& sf = str_fun->operator()(xbj, q2);
      cout << q2 << "|" << xbj << "|" << sf << endl;
      gr_f2.SetPoint(idx, lxbj, lq2, str_fun->F2(xbj, q2));
      //gr_f2.SetPoint( idx, xbj, q2, sf.F2 );
      gr_fl.SetPoint(idx, lxbj, lq2, str_fun->FL(xbj, q2));
      //gr_fl.SetPoint( idx, xbj, q2, sf.FL );
      idx++;
    }
  }

  //auto f_w2 = new TF1( "f_w2", "log10((10**x)/(1.-(10**x))*([0]-[1]**2))", lxbj_min, lxbj_max );
  /*auto f_w2 = new TF1( "f_w2", "log10((10**y)/(1.-(10**y))*([0]-[1]**2))", lxbj_min, lxbj_max );
  f_w2->SetParameter( 1, ParticleProperties::mass( PDG::Proton ) );
  f_w2->SetLineColor( kGreen+1 );
  f_w2->SetRange( lxbj_min, lxbj_max );*/
  /*  auto f_q2 = new TF1( "f_q2", "[0]", lxbj_min, lxbj_max );
  f_q2->SetParameter( 0, log10( 9. ) );
//  f_q2->SetParameter( 0, 0.1 );
  f_q2->SetLineColor( kGray+1 );
  f_q2->SetRange( lxbj_min, lxbj_max );
  cout << f_q2->Eval( lxbj_min ) << "|" << f_q2->Eval( lxbj_max ) <<endl;*/
  //  TF1 f_q2( "f_q2", "log10( 9. )", lxbj_min, lxbj_max );
  gStyle->SetPalette(kColorPrintableOnGrey);
  //gStyle->SetPalette( kSunset );
  //gStyle->SetPalette( 1 );
  //gStyle->SetPalette( 55 );
  gStyle->SetPalette(kLightTemperature);
  //gStyle->SetPalette( kGreyScale );
  //gStyle->SetPalette( kDeepSea );
  //gStyle->SetPalette( kTemperatureMap );
  //gStyle->SetPalette( kValentine );
  //gStyle->SetPalette( kGistEarth );
  //gStyle->SetPalette( kViridis );
  TColor::InvertPalette();
  //gStyle->SetNdivisions( 99, "Z" );
  std::ostringstream os;
  //os << type;
  os << StructureFunctionsFactory::get().describe(type);
  const char* plot_type =
      //    "cont1";
      //    "tri2";
      //    "pcol";
      //    "cont4z";
      //    "colz";
      "surf3";
  //    "scat";
  //    "lego2 0";
  {
    //Canvas c( "test_f2", Form( "%s F_{2}^{p} structure functions", os.str().c_str() ) );
    Canvas c(Form("%s_f2", basename), Form("%s F_{2}^{p} structure function", os.str().c_str()));
    //c.SetLogx();
    //c.SetLogy();
    gr_f2.Draw(plot_type);
    gr_f2.GetHistogram()->SetTitle(";log_{10}(x_{Bj});log_{10}(Q^{2}/GeV^{2});F_{2}^{p}(x_{Bj},Q^{2})");
    //gr_f2.GetHistogram()->SetContour( 90 );
    //gr_f2.SetMaximum( 6. );
    //c.SetLogz();
    /*f_w2->SetParameter( 0, 3. );
    f_w2->DrawCopy( "same" );
    f_w2->SetParameter( 0, 4. );
    f_w2->DrawCopy( "same" );
    f_q2.DrawCopy( "same");*/
    c.Prettify(gr_f2.GetHistogram());
    gr_f2.GetXaxis()->SetTitleOffset(1.4);
    gr_f2.GetYaxis()->SetTitleOffset(1.8);
    gr_f2.GetZaxis()->SetTitleOffset(1.2);
    //gr_f2.GetXaxis()->SetRangeUser( xbj_min, xbj_max );
    //gr_f2.GetYaxis()->SetRangeUser( q2_min, q2_max );
    //gPad->SetPhi( M_PI/3 );
    c.Save("png");
  }
  {
    //Canvas c( "test", "LUXlike F_{2}^{p}(x_{Bj},Q^{2})/Q^{2}" );
    //Canvas c( "test_fl", Form( "%s F_{L}^{p} structure functions", os.str().c_str() ) );
    Canvas c(Form("%s_fl", basename), Form("%s F_{L}^{p} structure function", os.str().c_str()));
    gr_fl.Draw(plot_type);
    //    gr_fl.Draw( "lego2 0" );
    gr_fl.GetHistogram()->SetTitle(";log_{10}(x_{Bj});log_{10}(Q^{2}/GeV^{2});F_{L}^{p}(x_{Bj},Q^{2})");
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
