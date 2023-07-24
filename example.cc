
#include "TApplication.h"

#include "iostream"
#include "iomanip"

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TRandom.h"

#include <gsl_sf_hyperg.h>

using namespace std;

Double_t _gaussexp( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];

  Double_t Norm = par[0];
  Double_t Q = par[1];
  Double_t s = par[2];
  Double_t alpha = par[3];
  Double_t w = par[4];
  
  Double_t arg = 0.0; 
  
  if ( s!=0.0 ) arg = ( xx - Q )/s;    
  else cout << "Error: The code tries to divide by zero." << endl;

  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  
  Double_t result = 0.0;
  if ( xx>=0.0 ) result = w*alpha*TMath::Exp( -xx*alpha ) + ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s*gn )*TMath::Exp( -0.5*arg*arg );

  result *= Norm;
  
  return result;
  
}

Float_t wbin;
Int_t nlim = 10;

Double_t _func( Double_t *x, Double_t *par )
{
  Double_t result = 0.0; 

  Double_t xx = x[0];
   
  Double_t Norm = par[0];
  
  Double_t Q0 = par[1];
  Double_t s0 = par[2];
  
  Double_t mu = par[3];
    
  Double_t Q = par[4];
  Double_t s = par[5];
  
  Double_t alpha = par[6];
  Double_t w = par[7];

  /* ... */

  Double_t arg = 0.0; 
  if ( s0!=0.0 ) arg = ( xx - Q0 )/s0;    
  else cout << "Error: The code tries to divide by zero." << endl;
    
  Double_t SR0 = 1.0/( sqrt( 2.0*TMath::Pi() )*s0 )*TMath::Exp( -0.5*arg*arg );
  SR0 *= TMath::Poisson( 0, mu );
  result += SR0; // 0
  
  
  Double_t omega0 = ( xx - Q0 - alpha*pow( s0, 2.0 ) )/sqrt(2.0)/s0;
  Double_t SR1 = w*alpha/2.0*TMath::Exp( -alpha*( xx-Q0 )+pow( alpha*s0, 2.0 )/2.0 )*TMath::Erfc( -omega0 );

  Double_t Q1 = Q0+Q;
  Double_t s12 = pow( s0, 2.0 )+pow( s, 2.0 );
  Double_t s1 = sqrt( s12 );
  
  Double_t arg1 = 0.0; 
  if ( s1!=0.0 ) arg1 = ( xx - Q1 )/s1;    
  else cout << "Error: The code tries to divide by zero." << endl;

  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  Double_t A = ( Q0-xx )*pow( s, 2.0 ) - Q*pow( s0, 2.0 ); 
  Double_t B = sqrt( 2.0 )*s0*s*s1;
  SR1 += ( 1.0-w )/2.0/gn/( sqrt( 2.0*TMath::Pi() )*s1 )*TMath::Exp( -0.5*arg1*arg1 )*TMath::Erfc( A/B );
  SR1 *= TMath::Poisson( 1, mu );
  result += SR1; // 1
   

  Double_t k = s/gn/sqrt( 2.0*TMath::Pi() )*TMath::Exp( -pow( Q, 2.0 )/( 2.0*pow( s, 2.0 ) ) );
  Double_t Qg = Q + k;
  Double_t sg2 = pow( s, 2.0 ) - ( Q+k )*k;
  
  for ( Int_t n = 2; n<nlim; n++ )
    {
      Double_t SRn = 0.0;
            
      Double_t Qn = Q0 + 1.0*n*Qg;
      Double_t sn2 = pow( s0,2.0 )+ 1.0*n*sg2;
      Double_t sn = sqrt( sn2 );
      
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx - Qn )/sn;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t gnB = 1.0/( sqrt( 2.0*TMath::Pi() )*sn )*TMath::Exp( -0.5*argn*argn );
      SRn += pow( 1.0-w, n )*gnB;
      
      
      for ( Int_t m=1; m<=n; m++ )
	{
	  Double_t Qmn = Q0 + 1.0*(n-m)*Qg;
	  Double_t smn2 = pow( s0, 2.0 )+1.0*(n-m)*sg2;
	  Double_t smn = sqrt( smn2 );
	  
	  Double_t cmn = alpha*pow( alpha*smn*sqrt( 2.0 ), m-1.0 )/TMath::Factorial( m-1.0 )/2.0/sqrt( TMath::Pi() );
	  
	  Double_t psi = ( xx-Qmn )/sqrt(2.0)/smn;
	  Double_t psi2 = pow( psi, 2.0 );
	  Double_t omega = ( xx-Qmn-alpha*pow( smn, 2.0 ) )/sqrt(2.0)/smn;
	  Double_t omega2 = pow( omega, 2.0 );
	  
	  Double_t A1m = 1.0*m/2.0;
	  Double_t A2m = (0.0+m+1.0)/2.0;
	  
	  Double_t Imn=0.0;
	  Double_t hi_limit=25.0;
	  
	  if ( omega>=hi_limit )
	    {
	      Imn = 2.0*sqrt( TMath::Pi() )*TMath::Exp( omega2-psi2 + ( m-1.0 )*TMath::Log( omega ) );
	      
	    }
	  else if ( omega<hi_limit && omega>=0.0  )
	    {
	      Double_t t1 = TMath::Gamma( A1m )*gsl_sf_hyperg_1F1( 1.0/2.0-A1m, 1.0/2.0, -omega2 );
	      Double_t t2 = 2.0*omega*TMath::Gamma( A2m )*gsl_sf_hyperg_1F1( 3.0/2.0-A2m, 3.0/2.0, -omega2 );
	      Imn = ( t1+t2 )*TMath::Exp( omega2-psi2 );
	      
	    }
	  else if ( omega<0.0 )
	    {
	      Double_t t3 = TMath::Gamma( A1m )*TMath::Gamma( A2m )/sqrt( TMath::Pi() );
	      Imn = t3*gsl_sf_hyperg_U( A1m, 1.0/2.0, omega2 )*TMath::Exp( -psi2 );
	      
	    }
	  
	  
	  Double_t hmnB = cmn*Imn;
	  Double_t binom = TMath::Factorial( n )/TMath::Factorial( m )/TMath::Factorial( n-m );
	  SRn += binom*pow( w, m )*pow( 1.0-w, n-m )*hmnB;
	  
	}
           
      SRn *= TMath::Poisson( n, mu );
      result += SRn; // n = 2-nlim
      
    } 
  
  
  Double_t Qs = w/alpha + (1.0-w)*Qg;
  Double_t ss2 = w/pow( alpha, 2.0 ) + (1-w)*sg2 + w*(1.0-w)*pow( Qg-1.0/alpha, 2.0 );
      
  for ( Int_t n = nlim; n<65; n++ )
    {
      Double_t Qn = Q0 + 1.0*n*Qs;
      Double_t sn2 = pow( s0, 2.0 ) + 1.0*n*ss2;
      Double_t sn = sqrt( sn2 );
      
      Double_t argn = 0.0; 
      if ( sn!=0.0 ) argn = ( xx - Qn )/sn;    
      else cout << "Error: The code tries to divide by zero." << endl;
      Double_t SRn = 1.0/( sqrt( 2.0*TMath::Pi() )*sn )*TMath::Exp( -0.5*argn*argn );

      SRn *= TMath::Poisson( n, mu );
      result += SRn; // n >= nlim
      
    }
      
  

  /* ... */
    
  result *= Norm*wbin;
  
  return result;
  
}

Int_t example()
{  
  time_t start;  
  time( &start );

  cout << "" << endl;
  cout << " the macro now begins ... ( example.cc )" << endl;
  cout << "" << endl;
  
  gROOT->Reset();
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetCanvasColor( 0 );
  gStyle->SetPadColor( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetFrameBorderMode( 0 );
  gStyle->SetTitleColor( 0 );   
  gStyle->SetTitleFillColor( 0 );  
  gStyle->SetTitleBorderSize( 0 );
  gStyle->SetTitleX( 0.10 );
  gStyle->SetTitleY( 0.98 );
  gStyle->SetTitleFont( 22, "" );
  gStyle->SetTitleSize( 0.055, ""  );
  gStyle->SetStatColor( 0 );
  gStyle->SetStatFont( 22 );
  gStyle->SetStatBorderSize( 1 );
  gStyle->SetStatX( 0.90 );
  gStyle->SetStatY( 0.90 );
  gStyle->SetStatFontSize( 0.04 );
  gStyle->SetOptStat( 1110 );
  gStyle->SetTitleFont( 22, "XYZ"  );
  gStyle->SetTitleSize( 0.05, "XYZ"  );
  gStyle->SetTitleColor( kBlack, "XYZ"  );
  gStyle->SetTitleAlign(13);
  gStyle->SetLabelFont( 22, "XYZ"  );
  gStyle->SetLabelSize( 0.04, "XYZ"  );
  gStyle->SetOptStat( 1110 ); 
  gStyle->SetOptFit( 0 ); 
  
  gStyle->SetPalette( 1 ); 
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
  gStyle->SetNumberContours( NCont );

  gROOT->ForceStyle();
    
  
  // Input <--------------------- !!!

  Float_t xmin = -200.0;
  Float_t xmax = 1300.0;
  Int_t nbins = 500;
  
  Int_t ntot = 2.5e+6;
  
  Float_t Q0 = 0.0;
  Float_t s0 = 21.0;
    
  Float_t mu = 1.5;
  
  Float_t Q = 120.0;
  Float_t s = 35.0;
  Float_t alpha = 1.0/30.0;
  Float_t w = 0.2;
  
    
  // <--------------------------- !!!
  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  c1->SetLogy();
  
  TF1 *spe = new TF1( "spe" , _gaussexp, 0.0, xmax, 5 );
  spe->SetParNames( "Norm", "Q", "#sigma", "#alpha", "w" );
  spe->SetNpx( 2000 );
  spe->SetLineColor( kRed );
  spe->SetLineWidth( 2 );
  
  spe->SetParameters( 1.0, Q, s, alpha, w );
   
  
  TH1D *h1 = new TH1D( "h1", "; Charge [AU]; Entries", nbins, xmin, xmax );
  gRandom->SetSeed( 0 );
  
  for ( Int_t i=0; i<ntot; i++ )
    {
      Double_t QQ = gRandom->Gaus( Q0, s0 );
      
      Int_t n = gRandom->Poisson( mu );
      for ( Int_t j=0; j<n; j++ ) QQ += spe->GetRandom();

      h1->Fill( QQ );
      
    }

  h1->SetMarkerStyle( 20 );
  h1->SetMarkerSize( 0.65 );
  h1->SetLineColor( kBlack );
  h1->SetLineWidth( 2.0 );
  h1->SetMarkerColor( kBlack );

  h1->Draw( "PEZ" );
  

  TF1 *ana = new TF1( "ana" , _func, xmin, xmax, 8 );
  ana->SetParNames( "Norm", "Q_{0}", "#sigma_{0}", "#mu", "Q", "#sigma", "#alpha", "w" );
  ana->SetNpx( 2000 );
  ana->SetLineColor( kAzure+6 );
  ana->SetLineWidth( 2 );
  
  ana->SetParameters( ntot, Q0, s0, mu, Q, s, alpha, w );
  wbin = h1->GetBinWidth( 1 );

  h1->Fit( "ana", "", "", xmin, xmax );
    
  c1->Update();
  c1->WaitPrimitive();
  
  cout << "" << endl;
  cout << "" << endl;
  
  time_t end;      
  time( &end );
      
  Int_t dura = difftime( end, start );      
  Int_t min = dura / 60; Int_t sec = dura % 60;
    
  cout << " ---> "<< Form( "%02d:%02d", min, sec ) << endl;  
  cout << "" << endl;
  
  return 0;
  
}

int main()
{
  TApplication theApp( "App", 0, 0 );
    
  example();

}
