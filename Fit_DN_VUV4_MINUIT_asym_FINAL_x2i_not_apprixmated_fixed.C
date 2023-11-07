
//-------------------------------------------------------------------
//
//  Hamamatsu VUV4
//
//  This macro fit the Dark Noise
//
//
//  I assume Pe and Ph are temperature independent.
//  I simply rescale considering the shifting of the breakdwon voltage with temperature.
//
//
//  The PDE has been measured at -40 C so transaltion is 0 for the probability
//  0.16 depends from the PDE zero
//
//  giacomo@triumf.ca
//
//-------------------------------------------------------------------



#include "../../Colors.c"
#include "../../Markers.c"


//-------------------------------------------------------------------
// Constants
//-------------------------------------------------------------------

#define q 1.6022e-19 //[C] Elementary charge

#define kb 1.380e-23 //[J/K]  Boltzmann Constant

#define mel 9.109383e-31 //[Kg] Electron Mass

#define meff_el 2.27735e-31//[Kg]  Electron Mass meff_el=0.25*mel

#define hcut 1.05457e-34  //[J*s]   Cut Planck Constant

#define h 6.62607e-34  //[J*s]   Planck Constant

#define kbeV 8.61673324e-5 //[eV/K]  Boltzmann Constant


#define NFILE 5

#define NUMBERFILE_TO_FIT 5

string filePath = "/Users/giacomo/Documents/TRIUMF/GIT/data/HAMAMATSU_SiPM/VUV4_DATA_PUBLISHED/DN_AP_JAN2018/";

string filename[NFILE] = {
                          "minus_40",
                          "minus_60",
                          "minus_80",
                          "minus_95",
                          "minus_110"
                         };

using namespace std;


//-------------------------------------------------------------------
// TGraph with data
//
//  0: -40
//  1: -60
//  2: -80
//  3: -95
//  4: -110
//
//-------------------------------------------------------------------

TGraphErrors *G1[NUMBERFILE_TO_FIT];

//always 7 point so no need to have one for each graph
const int nPoints0 = 7;



//-------------------------------------------------------------------
//
// Fit functions
//
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// Equation to fit the SRH rate
//-------------------------------------------------------------------

Double_t fit_exp(Double_t *x,Double_t *par) {
    
    //Parameters
    
    double tau0=par[0];
    double Eg0=par[1];
    double traps_energy=par[2];
    double exponent=par[3];
   
    
    //Temperature
    Double_t one_over_T=x[0]; //1/T
    
    //Find_T
    Double_t T=1.0/x[0];
    
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;
    mp_star=mp_star*mel;
    
    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    Double_t alpha=4.73e-4;
    Double_t beta=636;

    Double_t Eg=Eg0-((alpha*T*T)/(T+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    Double_t carrier_lifetime=tau0*TMath::Power((T/300),exponent);//s
    //Double_t carrier_lifetime=tau0;//s
    
    Double_t csrh_1=ni/(2.0*carrier_lifetime*cosh(traps_energy/(kbeV*T)));//cm^-3 s^-1
   
    
    //std::cout<<" "<<traps_energy<<" "<<kbeV<<" "<<" "<<T<<" "<<kbeV*T<<" "<<traps_energy/(kbeV*T)<<" "<<csrh_1<<std::endl;
   
    
    Double_t csrh_total=(csrh_1)*1e6;//m^-3 s^-1
    
     return csrh_total;
    
}




//-------------------------------------------------------------------
// Equation to fit the SRH rate
//-------------------------------------------------------------------

Double_t theory_shockley(Double_t *x,Double_t *par) {
    
    //Temperature
    Double_t one_over_T=x[0]; //1/T
    
    //Find_T
    Double_t Temperature=1.0/x[0];
    
    //eV/K
    Double_t alpha=4.73e-4;
    //K
    Double_t beta=636;
    //Nc and Nv cm^-3
    Double_t Nc=2.8e19;
    Double_t Nv=1.04e19;
    Double_t tau=1.0e-7;//s
   

    //https://ecee.colorado.edu/~bart/book/eband5.htm
    Double_t Eg=1.166-((alpha*Temperature*Temperature)/(Temperature+beta));
    
    //This is in 1/(cm3)
    Double_t nie_cm3=TMath::Sqrt(Nc*Nv*TMath::Exp(-(Eg/(kbeV*Temperature))));
    
    Double_t nie_m3=nie_cm3/1.0e-6;
   
    
    //csrh_expected
    Double_t csrh_expected=nie_m3/(2.0*tau);
    
    
     return csrh_expected;
    
}


//-------------------------------------------------------------------
// Equation for fitting -40 C
//-------------------------------------------------------------------

Double_t Fit_40(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias --> Common
    Double_t fDN_1=par[1]; //Junction fraction  --> Common
    Double_t csrh=par[2]; //Schokley Read Hall  --> Not Common
    Double_t T=par[3]; //Temperature K          --> Not Common
    Double_t p=par[4]; //junction type          --> Common
    Double_t Vint=par[5]; //built-in            --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper     --> Common
    Double_t tr=par[7]; //translation for shift in the breakdwon due to the temperature   --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
        
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left.
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //std::cout<<" Gamma_trap "<<Gamma_trap<<std::endl;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1+term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}



//-------------------------------------------------------------------
// Equation for fitting -60 C
//-------------------------------------------------------------------

Double_t Fit_60(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias  --> Common
    Double_t fDN_1=par[1]; //Junction fraction   --> Common
    Double_t csrh=par[8]; //Schokley Read Hall   --> Not Common
    Double_t T=par[9]; //Temperature K           --> Not Common
    Double_t p=par[4]; //junction type           --> Common
    Double_t Vint=par[5]; //built-in             --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper   --> Common
    Double_t tr=par[10]; //translation for shift in the breakdwon due to the temperature --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1+term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}


//-------------------------------------------------------------------
// Equation for fitting -80 C
//-------------------------------------------------------------------

Double_t Fit_80(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias   --> Common
    Double_t fDN_1=par[1]; //Junction fraction    --> Common
    Double_t csrh=par[11]; //Schokley Read Hall   --> Not Common
    Double_t T=par[12]; //Temperature K           --> Not Common
    Double_t p=par[4]; //junction type            --> Common
    Double_t Vint=par[5]; //built-in              --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper    --> Common
    Double_t tr=par[13]; //translation for shift in the breakdwon due to the temperature --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    //the 0.16 may come from error breakdown. What about PDE ?
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1+term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}



//-------------------------------------------------------------------
// Equation for fitting -95 C
//-------------------------------------------------------------------

Double_t Fit_95(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias  --> Common
    Double_t fDN_1=par[1]; //Junction fraction   --> Common
    Double_t csrh=par[14]; //Schokley Read Hall  --> Not Common
    Double_t T=par[15]; //Temperature K          --> Not Common
    Double_t p=par[4]; //junction type           --> Common
    Double_t Vint=par[5]; //built-in             --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper    --> Common
    Double_t tr=par[16]; //translation for shift in the breakdwon due to the temperature --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    //the 0.16 may come from error breakdown. What about PDE ?
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1+term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
}



//-------------------------------------------------------------------
// Equation for fitting -110 C
//-------------------------------------------------------------------

Double_t Fit_110(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias   --> Common
    Double_t fDN_1=par[1]; //Junction fraction    --> Common
    Double_t csrh=par[17]; //Schokley Read Hall   --> Not Common
    Double_t T=par[18]; //Temperature K           --> Not Common
    Double_t p=par[4]; //junction type            --> Common
    Double_t Vint=par[5]; //built-in              --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper     --> Common
    Double_t tr=par[19]; //translation for shift in the breakdwon due to the temperature  --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    //the 0.16 may come from error breakdown. What about PDE ?
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    
    //std::cout<<" Gamma_trap "<<Gamma_trap<<std::endl;
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power(numerator/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1+term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
}

//-------------------------------------------------------------------
//
// Function to show the different contribution. The trick is to keep
// the same parameter naming
//
//-------------------------------------------------------------------



//-------------------------------------------------------------------
// Equation for fitting -40 C enhanced
//-------------------------------------------------------------------

Double_t Fit_40_enhanced(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias --> Common
    Double_t fDN_1=par[1]; //Junction fraction  --> Common
    Double_t csrh=par[2]; //Schokley Read Hall  --> Not Common
    Double_t T=par[3]; //Temperature K          --> Not Common
    Double_t p=par[4]; //junction type          --> Common
    Double_t Vint=par[5]; //built-in            --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper     --> Common
    Double_t tr=par[7]; //translation for shift in the breakdwon due to the temperature   --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    Double_t V_tr=x[0]+tr+0.16;  //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    //std::cout<<" -Gamma_trap "<<-Gamma_trap<<std::endl;
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}


//-------------------------------------------------------------------
// Equation for fitting -40 C standard
//-------------------------------------------------------------------

Double_t Fit_40_standard(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias --> Common
    Double_t fDN_1=par[1]; //Junction fraction  --> Common
    Double_t csrh=par[2]; //Schokley Read Hall  --> Not Common
    Double_t T=par[3]; //Temperature K          --> Not Common
    Double_t p=par[4]; //junction type          --> Common
    Double_t Vint=par[5]; //built-in            --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper     --> Common
    Double_t tr=par[7]; //translation for shift in the breakdwon due to the temperature   --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    Double_t V_tr=x[0]+tr+0.16;  //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}



//-------------------------------------------------------------------
// Equation for fitting -110 C enhanced
//-------------------------------------------------------------------

Double_t Fit_110_enhanced(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias   --> Common
    Double_t fDN_1=par[1]; //Junction fraction    --> Common
    Double_t csrh=par[17]; //Schokley Read Hall   --> Not Common
    Double_t T=par[18]; //Temperature K           --> Not Common
    Double_t p=par[4]; //junction type            --> Common
    Double_t Vint=par[5]; //built-in              --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper     --> Common
    Double_t tr=par[19]; //translation for shift in the breakdwon due to the temperature  --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    //the 0.16 may come from error breakdown. What about PDE ?
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_2+term_4);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}

//-------------------------------------------------------------------
// Equation for fitting -110 C enhanced
//-------------------------------------------------------------------


Double_t Fit_110_standard(Double_t *x,Double_t *par) {

    
    //-------------------------------------------------------------------
    // 1+2 Parameters floated and fixed
    //-------------------------------------------------------------------
    
    Double_t W0=par[0]; //Junction at zero bias   --> Common
    Double_t fDN_1=par[1]; //Junction fraction    --> Common
    Double_t csrh=par[17]; //Schokley Read Hall   --> Not Common
    Double_t T=par[18]; //Temperature K           --> Not Common
    Double_t p=par[4]; //junction type            --> Common
    Double_t Vint=par[5]; //built-in              --> Common
    Double_t k=par[6]; //effective ratio ionization coefficients ---> Avalanche paper     --> Common
    Double_t tr=par[19]; //translation for shift in the breakdwon due to the temperature  --> Not Common
    Double_t cbbt=par[20]; //cbbt
    
    Double_t fDN=fDN_1; //Junction fraction
    
    //-------------------------------------------------------------------
    // 3. x
    //-------------------------------------------------------------------
    
    Double_t V=x[0]; //Voltage not translated
    
    //the 0.16 may come from error breakdown. What about PDE ?
    Double_t V_tr=x[0]+tr+0.16; //Translation of the voltage to the left
    
    //-------------------------------------------------------------------
    // 4. Extrapolations
    //
    // The idea is to translate the probabilities to the left in order to accomodate
    // for shifting of the breakdown voltage
    //
    //-------------------------------------------------------------------
    
    
    Double_t Fm0=(Vint/((1-p)*W0));
    Double_t Fm=Fm0*TMath::Power(1+(V/Vint),1-p);
    Double_t W=W0*TMath::Power(1+V/Vint,p);
    
    //-------------------------------------------------------------------
    // Pe
    //-------------------------------------------------------------------
    
    Double_t ke=2.51354e+10; // Avalanche paper
    Double_t ke2=1.93824e+02; // Avalanche paper
    Double_t gamma_Pe=2.0; // Avalanche paper
    Double_t Exp_Vtot=0.5; // Avalanche paper
    Double_t gamma=ke*V_tr*TMath::Exp(-ke2/pow(V_tr,Exp_Vtot));
    Double_t Pe = (1.0-pow(gamma,-gamma_Pe));

    //-------------------------------------------------------------------
    // Ph --> From avalanche paper
    //-------------------------------------------------------------------
    
    Double_t Ph=(1.0-pow(1.0-Pe,k));
    
    
    //-------------------------------------------------------------------
    // Fit Function
    //-------------------------------------------------------------------
    
    
    Double_t epsilon0=8.854e-12;//C^2*N^-1*m^-2
    Double_t epsilonr=11.68;//#
    Double_t epsilon=epsilon0*epsilonr;
        
    Double_t mn_star=(1.045)+(0.00045*T);
    Double_t mp_star=(0.523)+(0.0014*T)-(1.48e-6)*T*T;
    
    mn_star=mn_star*mel;//Kg
    mp_star=mp_star*mel;//Kg
    

    Double_t Nc=2*TMath::Power(((2*TMath::Pi()*mn_star*kb*T)/(h*h)),(3.0/2.0));
    Double_t Nv=2*TMath::Power(((2*TMath::Pi()*mp_star*kb*T)/(h*h)),(3.0/2.0));
    
    //std::cout<<"Nc [m^-3]"<<Nc<<std::endl;//m^-3
    //std::cout<<"Nv [m^-3]"<<Nv<<std::endl;//m^-3
    
    Double_t Nc_cm=Nc*1e-6;//cm^-3
    Double_t Nv_cm=Nv*1e-6;//cm^-3
    
    //std::cout<<"Nc_cm [cm^-3]"<<Nc_cm<<std::endl;//m^-3
    //std::cout<<"Nv_cm [cm^-3]"<<Nv_cm<<std::endl;//m^-3

    
    Double_t alpha=4.73e-4;
    Double_t beta=636;
    
    Double_t Eg=1.16-((alpha*T*T)/(T+beta));//eV
    Double_t Eg_300=1.16-((alpha*300*300)/(300+beta));//eV
    
    Double_t ni=TMath::Sqrt(Nc_cm*Nv_cm)*TMath::Exp(-Eg/(2.0*kbeV*T));//cm^-3
    
    ni=ni*1e6;//m^-3
    
    Double_t Gamma_trap=TMath::Sqrt((kb*T*(W0*W0/Vint)/q)*TMath::Log((2*epsilon)/(q*((W0*W0)/(Vint))*ni)));
    
    //Gamma_trap=0;
    
    
    Double_t FGamma=TMath::Sqrt(24.0*meff_el*TMath::Power((kb*T),3.0))/(q*hcut);
    
    Double_t term_1=((W0*TMath::Sqrt(V))/TMath::Sqrt(Vint)-Gamma_trap);
    
    Double_t term_2=-((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm))*TMath::Exp(TMath::Power((Gamma_trap/((1-fDN)*W)),2.0));
    
    Double_t term_3=((TMath::Sqrt(3*TMath::Pi())*(1-fDN)*W*FGamma)/TMath::Abs(Fm));
    
    
    Double_t numerator=Fm*(1.0-(((1-fDN)*W-term_1-Gamma_trap)/((1-fDN)*W)));
    
    Double_t term_4=term_3*(TMath::Exp(TMath::Power((numerator)/FGamma,2.0)));
    
    Double_t Rdn=csrh*Ph*(term_1);
    
    //-------------------------------------------------------------------
    // BBT Term
    //-------------------------------------------------------------------
  
    

    Double_t F0=1.9e9 ;//V/m
    Double_t Rbbt=cbbt*V*TMath::Power(Fm,1.5)*TMath::Exp(-(F0/Fm)*(TMath::Power(Eg/Eg_300,1.5)));
    
    Rbbt=Rbbt*Pe;
    
    //-------------------------------------------------------------------
    // Final Fitting equation
    //-------------------------------------------------------------------
    
    
    Double_t DeltaE_eV=1.12/2.0;//eV
    Double_t condition=FGamma*TMath::Sqrt(DeltaE_eV/(3.0*kbeV*T));
    
    
    //Both holds
    if(Fm<condition){
       
        Rdn=Rbbt+Rdn;
        
    //only BBT
    }else{
        
        Rdn=Rbbt;
        
    }
    
    
    return Rdn;
    
}


//-------------------------------------------------------------------
//
// Chi Square
//
//-------------------------------------------------------------------


//-------------------------------------------------------------------
// Compute Chi Square -40 C
//-------------------------------------------------------------------
Double_t fcn0(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    const Int_t nbins = nPoints0;
    Int_t i;
    
    
    //calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
        Double_t delta  = ((G1[0]->GetY()[i])-Fit_40(&(G1[0]->GetX()[i]),par))/(G1[0]->GetErrorY(i));
        
        chisq += delta*delta;
        
        //cout<<chisq<<endl;
        
    }
    
    //cout << "FCN1 CHI2 " << chisq<<endl;
    
    return chisq;
    
}



//-------------------------------------------------------------------
// Compute Chi Square -60 C
//-------------------------------------------------------------------
Double_t fcn1(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    const Int_t nbins = nPoints0;
    Int_t i;
    
    
    //calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
        Double_t delta  = ((G1[1]->GetY()[i])-Fit_60(&(G1[1]->GetX()[i]),par))/(G1[1]->GetErrorY(i));
        
        chisq += delta*delta;
        
        //cout<<chisq<<endl;
        
    }
    
    //cout << "FCN1 CHI2 " << chisq<<endl;
    
    return chisq;
    
}



//-------------------------------------------------------------------
// Compute Chi Square -80 C
//-------------------------------------------------------------------
Double_t fcn2(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    const Int_t nbins = nPoints0;
    Int_t i;
    
    
    //calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
        Double_t delta  = ((G1[2]->GetY()[i])-Fit_80(&(G1[2]->GetX()[i]),par))/(G1[2]->GetErrorY(i));
        
        chisq += delta*delta;
        
        //cout<<chisq<<endl;
        
    }
    
    //cout << "FCN1 CHI2 " << chisq<<endl;
    
    return chisq;
    
}



//-------------------------------------------------------------------
// Compute Chi Square -95 C
//-------------------------------------------------------------------
Double_t fcn3(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    const Int_t nbins = nPoints0;
    Int_t i;
    
    
    //calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
        Double_t delta  = ((G1[3]->GetY()[i])-Fit_95(&(G1[3]->GetX()[i]),par))/(G1[3]->GetErrorY(i));
        
        chisq += delta*delta;
        
        //cout<<chisq<<endl;
        
    }
    
    //cout << "FCN1 CHI2 " << chisq<<endl;
    
    return chisq;
    
}


//-------------------------------------------------------------------
// Compute Chi Square -110 C
//-------------------------------------------------------------------
Double_t fcn4(Int_t &npar, Double_t *gin, Double_t *par, Int_t iflag)
{
    const Int_t nbins = nPoints0;
    Int_t i;
    
    
    //calculate chisquare
    Double_t chisq = 0;
    for (i=0;i<  nbins; i++) {
        Double_t delta  = ((G1[4]->GetY()[i])-Fit_110(&(G1[4]->GetX()[i]),par))/(G1[4]->GetErrorY(i));
        
        chisq += delta*delta;
        
        //cout<<chisq<<endl;
        
    }
    
    //cout << "FCN1 CHI2 " << chisq<<endl;
    
    return chisq;
    
}


//-------------------------------------------------------------------
// Compute Chi Square Combined Function
//-------------------------------------------------------------------

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    
    /*
    for(int i = 0 ; i < 5; i++) {
        std::cout << "P" << i << ": " << par[i] << std::endl;
    }
     */
    f = fcn0(npar,gin,par, iflag) +fcn1(npar,gin,par, iflag) + fcn2(npar,gin, par, iflag)+ fcn3(npar,gin, par, iflag)+ fcn4(npar,gin, par, iflag);
    
    //f = fcn1(npar,gin,par, iflag);
    
    //f = fcn2(npar,gin, par, iflag);
}






//-------------------------------------------------------------------
// Fit Routine
//-------------------------------------------------------------------

void fit_DN_multi_asym(){
    
    
    int aWhich,l;
    
    TCanvas* C1 = new TCanvas("C1","Graph with Fit",50,50,700,500);
    
    
    C1->SetBorderMode(0);
    C1->SetFillColor(0);
    C1->SetLeftMargin(0.12);
    C1->SetRightMargin(0.04);
    C1->SetBottomMargin(0.1);
    C1->SetGrid();
    C1->SetLogy();
    
    Double_t *Voltage=new Double_t[2000];
    Double_t *Voltage_err=new Double_t[2000];
    Double_t *DN=new Double_t[2000];
    Double_t *DN_err=new Double_t[2000];
    
    ifstream infile[NUMBERFILE_TO_FIT];
    
    //Fit alal files
    int data[NUMBERFILE_TO_FIT]={0,1,2,3,4};
    
    //Decleare global legend
    TLegend *leg2 = new TLegend(0.13,0.6,0.28,0.89);
    leg2->SetHeader("Temperature:");
    leg2->SetTextSize(gStyle->GetLabelSize());
    
    
    
    for(l=0;l<NUMBERFILE_TO_FIT;l++){
        
      
        
        //cout<<"l: "<<l<<endl;
    
        aWhich=data[l];
        
        cout<<"aWhich :"<<aWhich<<endl;
        
        
        Double_t Voltage_import,OV,OV_error,DN_inport,DN_Err_inport,AP,AP_err;
        
        infile[l].open((filePath+filename[aWhich]+".txt").c_str());
        
        if(!infile[l].is_open()) {
            std::cout<<" Filename: "<<filePath+filename[aWhich]+".txt Not Found !"<<std::endl;
            return;
        }
        
        Int_t cnt =0;
        cout<<"Importing data from File ..."<<filename[aWhich]<<endl;
        
        //Skip first dummy line with file info
        string dummyLine;
        getline(infile[l],dummyLine);
        
        
        while(!infile[l].eof())
        {
        
                            
            infile[l]>>Voltage_import>>OV>>OV_error>>DN_inport>>DN_Err_inport>>AP>>AP_err;
            
            cout<<Voltage_import<<" "<<DN_inport<<" "<<DN_Err_inport<<" "<<endl;
            
            Voltage[cnt]=Voltage_import;
        
            //From Datasheet 50 V Range
            Voltage_err[cnt]=(0.0015*Voltage_import+40e-3);
            
            DN[cnt]=DN_inport*1e6;//go to Hz/m2
            
            DN_err[cnt]=DN_Err_inport*1e6; // Error wasa bit underestiamted. Consider an extra 2% (is 2% by definition)
            
            DN[cnt]=DN[cnt]/0.6;//account fill factor
            
            
            //DN[cnt]=DN_inport;//go to Hz/m2
            //DN_err[cnt]=DN_Err_inport;
            
            cnt++;
        }
        
        cout<<"Number of lines "<<cnt<<endl;
        infile[l].close();
        

        //Increase number of cicles to fit
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000000);
        
        G1[l]=new TGraphErrors(cnt,Voltage,DN,Voltage_err,DN_err);

        if(l==0)   leg2->AddEntry(G1[l],"233 [K]","lp");
        if(l==1)   leg2->AddEntry(G1[l],"213 [K]","lp");
        if(l==2)   leg2->AddEntry(G1[l],"193 [K]","lp");
        if(l==3)   leg2->AddEntry(G1[l],"178 [K]","lp");
        if(l==4)   leg2->AddEntry(G1[l],"163 [K]","lp");
        
        G1[l]->SetMarkerColor(Colors(aWhich));
        G1[l]->SetLineColor(Colors(aWhich));
        G1[l]->SetMarkerStyle(Markers(aWhich));
        G1[l]->SetMarkerSize(1.4);
        //G1[l]->SetFillColorAlpha(Colors(aWhich),0.5);
        //G1[l]->SetLineStyle(9);
        G1[l]->SetTitle(" ");
        
    }
    
    
    TMultiGraph *mg1 = new TMultiGraph();
    
    mg1->SetTitle("");
    mg1->Add(G1[0]);
    mg1->Add(G1[1]);
    mg1->Add(G1[2]);
    mg1->Add(G1[3]);
    mg1->Add(G1[4]);
    // draw this MultiGraph
    mg1->Draw("AP");
    leg2->Draw("same");
    
    // X-axis
    TAxis *axisx1=mg1->GetXaxis();
    axisx1->SetTitle("voltage [V]");
    axisx1->SetRangeUser(45,57);
    axisx1->SetLimits(45,57);
    axisx1->SetTitleOffset(1.2);
    //axisx1->SetRangeUser(30,41);
    axisx1->Draw();
    
    // Y-axis
    TAxis *axisy1=mg1->GetYaxis();
    axisy1->SetTitle("dark noise [Hz/m^{2}]");
    axisy1->SetTitleOffset(1.2);
    axisy1->Draw();

        
    
    //---------------------------------------------------------------------
    //     Blocks for fitting procedure
    //---------------------------------------------------------------------
    
    //TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter *mini = TVirtualFitter::Fitter(0,21); //  TVirtualFitter *Fitter(TObject *obj, Int_t maxpar = 25);
    mini->SetFCN(fcn);
    Double_t arglist[100];
    arglist[0] = 2000000000; // number of function calls
    arglist[1] = 0.001; // tolerance
    

    
    Double_t vstart[21] = {
        2.00941e-07,             //W0
        0,             //fDN
        3.70367e+12,          //CSRH_40
        233,         //T_40
        0.5,       //p
        0.7,       //Vint
        0.07,       //k
        0,       //TR_40
        7.70367e+10,       //CSRH_60
        213,       //T_60
        1.0148,       //TR_60
        3.70367e+09,       //CSRH_80
        193,       //T_80
        2.1481,       //TR_80
        1e+08,       //CSRH_95
        178,       //T_95
        2.794,       //TR_95
        3.70367e7,       //CSRH_110
        163,       //T_100
        3.4568,       //TR_100
        3.38254e+12
    };
    
    
    
    
    Double_t vstart_min[21]   = {
        1e-8, //W0
        0.0001,    //fDN
        1e5,  //CSRH_40
        0,    //T_40
        0,    //p
        0,    //Vint
        0,    //k
        vstart[7]-0.2,    //TR_40
        1e1,  //CSRH_60
        0,    //T_60
        vstart[10]-0.2,    //TR_60
        1e1,  //CSRH_80
        0,    //T_80
        vstart[13]-0.2,    //TR_80
        1e1,  //CSRH_95
        0,    //T_95
        vstart[16]-0.2,    //TR_95
        1e1,  //CSRH_110
        0,    //T_100
        vstart[19]-0.2, //TR_100
        1e10
    };
    
    
    
    Double_t vstart_max[21]   = {
        1e-5, //W0
        0.9,    //fDN
        1e15,  //CSRH_40
        1,    //T_40
        1,    //p
        1,    //Vint
        1,    //k
        vstart[7]+0.2,    //TR_40
        1e15,  //CSRH_60
        1,    //T_60
        vstart[10]+0.2,    //TR_60
        1e15,  //CSRH_80
        1,    //T_80
        vstart[13]+0.2,    //TR_80
        1e15,  //CSRH_95
        1,    //T_95
        vstart[16]+0.2 ,    //TR_95
        1e15,  //CSRH_110
        1,    //T_100
        vstart[19]+0.2,      //TR_100
        1e17
    };
    
    
    Double_t step_2[21]={
        1e-6,
        1e-1,
        1e1,
        1,
        0.1,
        0.1,
        0.01,
        0.5,
        1e1,
        1,
        0.1,
        1e1,
        1,
        0.1,
        1e1,
        1,
        0.1,
        1e1,
        1,
        0.01,
        1
    };
    
    mini->ExecuteCommand("SET PRINT",arglist,21);
    mini->SetParameter(0,"W0", vstart[0],step_2[0],vstart_min[0],vstart_max[0]);
    mini->SetParameter(1,"fDN", vstart[1],step_2[1],vstart_min[1],vstart_max[1]);
    mini->SetParameter(2,"CSRH_40", vstart[2],step_2[2],vstart_min[2],vstart_max[2]);
    mini->SetParameter(3,"T_40", vstart[3],step_2[3],vstart_min[3],vstart_max[3]);
    mini->SetParameter(4,"p", vstart[4],step_2[4],vstart_min[4],vstart_max[4]);
    mini->SetParameter(5,"Vint", vstart[5],step_2[5],vstart_min[5],vstart_max[5]);
    mini->SetParameter(6,"k", vstart[6],step_2[6],vstart_min[6],vstart_max[6]);
    mini->SetParameter(7,"TR_40", vstart[7],step_2[7],vstart_min[7],vstart_max[7]);
    mini->SetParameter(8,"CSRH_60", vstart[8],step_2[8],vstart_min[8],vstart_max[8]);
    mini->SetParameter(9,"T_60", vstart[9],step_2[9],vstart_min[9],vstart_max[9]);
    mini->SetParameter(10,"TR_60", vstart[10],step_2[10],vstart_min[10],vstart_max[10]);
    mini->SetParameter(11,"CSRH_80", vstart[11],step_2[11],vstart_min[11],vstart_max[11]);
    mini->SetParameter(12,"T_80", vstart[12],step_2[12],vstart_min[12],vstart_max[12]);
    mini->SetParameter(13,"TR_80", vstart[13],step_2[13],vstart_min[13],vstart_max[13]);
    mini->SetParameter(14,"CSRH_95", vstart[14],step_2[14],vstart_min[14],vstart_max[14]);
    mini->SetParameter(15,"T_95", vstart[15],step_2[15],vstart_min[15],vstart_max[15]);
    mini->SetParameter(16,"TR_95", vstart[16],step_2[16],vstart_min[16],vstart_max[16]);
    mini->SetParameter(17,"CSRH_110", vstart[17],step_2[17],vstart_min[17],vstart_max[17]);
    mini->SetParameter(18,"T_110", vstart[18],step_2[18],vstart_min[18],vstart_max[18]);
    mini->SetParameter(19,"TR_110", vstart[19],step_2[19],vstart_min[19],vstart_max[19]);
    mini->SetParameter(20,"cbbt", vstart[20],step_2[20],vstart_min[20],vstart_max[20]);
    
    //mini->FixParameter(0);
    //mini->FixParameter(1);
    //mini->FixParameter(2);
    mini->FixParameter(3); //T_40
    mini->FixParameter(4); //p
    mini->FixParameter(5); //Vint
    mini->FixParameter(6); //k
    mini->FixParameter(7); //TR_40
    //mini->FixParameter(8);
    mini->FixParameter(9); //T_60
    mini->FixParameter(10); //TR_60
    //mini->FixParameter(11);
    mini->FixParameter(12); //T_80
    mini->FixParameter(13); //TR_80
    //mini->FixParameter(14);
    mini->FixParameter(15); //T_95
    mini->FixParameter(16); //TR_95
    //mini->FixParameter(17);
    mini->FixParameter(18); //T_110
    mini->FixParameter(19); //TR_110
    //mini->FixParameter(20); //CBBT
    
    
    mini->ExecuteCommand("MIGRAD",arglist,2);
    mini->ExecuteCommand("IMPROVE",arglist,2);
    
    Double_t we,al,bl, chiS;
    Char_t parName[64];
    
    
    Double_t currentPar[21] = {0};
    for (int i=0; i< 21;i++) {
      mini->GetParameter(i,parName,currentPar[i],we,al,bl);
    }
    
    
    
    
    
    
    
    TF1 *fun_0=new TF1("fun_0",Fit_40,45,57,21);
    
    fun_0->SetParameters(currentPar);
    fun_0->SetLineColor(Colors(0));
    fun_0->SetLineStyle(1);
    //fun_0->SetLineWidth(4);
    
    //     mg1->Add(new TGraph(fun_1));
    fun_0->Draw("same");
    
    
    
    TF1 *fun_1=new TF1("fun_1",Fit_60,45,57,21);
    
    fun_1->SetParameters(currentPar);
    fun_1->SetLineColor(Colors(1));
    fun_1->SetLineStyle(1);
    //fun_1->SetLineWidth(4);

    //     mg1->Add(new TGraph(fun_1));
    fun_1->Draw("same");
   
    TF1 *fun_2=new TF1("fun_2",Fit_80,45,57,21);
    fun_2->SetParameters(currentPar);
    fun_2->SetLineColor(Colors(2));
    fun_2->SetLineStyle(1);
    //fun_2->SetLineWidth(4);
    //     mg1->Add(new TGraph(fun_2));
    fun_2->Draw("same");
    
    TF1 *fun_3=new TF1("fun_3",Fit_95,45,57,21);
    fun_3->SetParameters(currentPar);
    fun_3->SetLineColor(Colors(3));
    fun_3->SetLineStyle(1);
    //fun_3->SetLineWidth(4);
    //     mg1->Add(new TGraph(fun_2));
    fun_3->Draw("same");
    
    
    TF1 *fun_4=new TF1("fun_4",Fit_110,45,57,21);
    fun_4->SetParameters(currentPar);
    fun_4->SetLineColor(Colors(4));
    fun_4->SetLineStyle(1);
    //fun_4->SetLineWidth(4);
    //     mg1->Add(new TGraph(fun_2));
    fun_4->Draw("same");

    
    
    //---------------------------------------------------------------------
    //    Extrapolate saturation
    //---------------------------------------------------------------------
    
    //GetParameter (Int_t ipar, char *name, Double_t &value, Double_t &verr, Double_t &vlow, Double_t &vhigh) const =0
    //Last two are just the parameter boundaries.
    
    TGraphErrors *gSRHall_rate[1];
    gSRHall_rate[0] = new TGraphErrors(); //SRH rate
    
    Double_t temp,temp_2;
    
    Double_t err,err_2;
    
    //Increase all the errors of a factor 11 of CSRH
    
    mini->GetParameter(2,parName,temp,err,al,bl); //CSRH -40 C
    mini->GetParameter(3,parName,temp_2,err_2,al,bl); //T -40 C
    
    gSRHall_rate[0]->SetPoint(0,1.0/temp_2,temp);
    gSRHall_rate[0]->SetPointError(0,0,err*11);
    
    
    mini->GetParameter(8,parName,temp,err,al,bl); //CSRH -60 C
    mini->GetParameter(9,parName,temp_2,err_2,al,bl); //T -60 C
    gSRHall_rate[0]->SetPoint(1,1.0/temp_2,temp);
    gSRHall_rate[0]->SetPointError(1,0,err*11);
    
    
    mini->GetParameter(11,parName,temp,err,al,bl); //CSRH -80 C
    mini->GetParameter(12,parName,temp_2,err_2,al,bl); //T -80 C
    gSRHall_rate[0]->SetPoint(2,1.0/temp_2,temp);
    gSRHall_rate[0]->SetPointError(2,0,err*11);
    
    mini->GetParameter(14,parName,temp,err,al,bl); //CSRH -95 C
    mini->GetParameter(15,parName,temp_2,err_2,al,bl); //T -95 C
    gSRHall_rate[0]->SetPoint(3,1.0/temp_2,temp);
    gSRHall_rate[0]->SetPointError(3,0,err*11);
    
    mini->GetParameter(17,parName,temp,err,al,bl); //CSRH -110 C
    mini->GetParameter(18,parName,temp_2,err_2,al,bl); //T -110 C
    gSRHall_rate[0]->SetPoint(4,1.0/temp_2,temp);
    gSRHall_rate[0]->SetPointError(4,0,err*11);
    
    
    TCanvas *CDplot2 = new TCanvas("c2","c2",1000,600);
    
    CDplot2->SetGrid();
    CDplot2->SetBorderMode(0);
    CDplot2->SetFillColor(0);
    CDplot2->SetLeftMargin(0.12);
    CDplot2->SetRightMargin(0.04);
    CDplot2->SetBottomMargin(0.1);
    CDplot2->Draw();

    gSRHall_rate[0]->SetName("SRH_data");
    
    gSRHall_rate[0]->GetXaxis()->SetTitle("1/T [1/K]");
    gSRHall_rate[0]->GetXaxis()->SetLimits(4e-3,6.5e-3);
    

    gSRHall_rate[0]->SetMarkerColor(Colors(0));
    gSRHall_rate[0]->SetLineColor(Colors(0));

    
    gSRHall_rate[0]->SetMarkerStyle(Markers(0));
    gSRHall_rate[0]->SetMarkerSize(1.2);
    
    gSRHall_rate[0]->GetYaxis()->SetTitle("SRH [Hz/m^{3}]");
    gSRHall_rate[0]->GetYaxis()->SetRangeUser(1e5,1e15);
    gSRHall_rate[0]->GetYaxis()->SetTitleOffset(1.2);
    gSRHall_rate[0]->GetXaxis()->SetTitleOffset(1.2);
    gSRHall_rate[0]->Draw("ap");
    CDplot2->SetLogy();
    
    //---------------------------------------------------------------------
    //    Fit saturation
    //---------------------------------------------------------------------
    
    
    
    TF1* fexp = new TF1("fit_exp",fit_exp,4.2e-3,6.4e-3,4);
         
    fexp->SetParName(0,"Lifetime 300K");
    fexp->SetParName(1,"Eg");
    fexp->SetParName(2,"traps_energy");
    fexp->SetParName(3,"exponents");
         
         
    //Sat
    fexp->FixParameter(0,2.5e-3);
    fexp->SetParLimits(0,1e-4,2.5e-3);
         
    //Eg
    fexp->FixParameter(1,1.166);
    
    //traps
    fexp->SetParameter(2,0);
    
    //traps exp
    fexp->FixParameter(3,+0.57);
    
         
    gSRHall_rate[0]->Fit("fit_exp","R");
    fexp->Draw("same,l");
    
    
    std::cout<<" FIT SRH ....."<<std::endl;
    
    std::cout<<" The Reduced Chi Square is "<<fexp->GetChisquare()/fexp->GetNDF()<<std::endl;
    
             
    
    //---------------------------------------------------------------------
    //    Now plot theoretical shockley value
    //---------------------------------------------------------------------
    
    /*
    TF1* fshock_expected = new TF1("fit_exp",theory_shockley,4.2e-3,6.4e-3,0);
        
    fshock_expected->SetLineColor(Colors(1));
    fshock_expected->SetLineStyle(5);
    fshock_expected->SetLineWidth(2);
    fshock_expected->Draw("same,l");
    */
    
    
    //---------------------------------------------------------------------
    //    Separate Contribution
    //---------------------------------------------------------------------
    
    //Decleare global legend
    TLegend* leg3 = new TLegend(0.8,0.4,0.9,0.6);
    
    
    
    TCanvas* C3 = new TCanvas("C3","Show contribution",50,50,700,500);
    
    
    C3->SetBorderMode(0);
    C3->SetFillColor(0);
    C3->SetLeftMargin(0.12);
    C3->SetRightMargin(0.04);
    C3->SetBottomMargin(0.1);
    C3->SetGrid();
    C3->SetLogy();
    C3->Draw();
    
    //-40 C
    G1[0]->GetYaxis()->SetRangeUser(1e-3,1e9);
    G1[0]->GetYaxis()->SetLimits(1e-3,1e9);
    G1[0]->GetXaxis()->SetRangeUser(44,57);
    G1[0]->GetXaxis()->SetLimits(44,57);
    G1[0]->GetYaxis()->SetTitle("Dark Noise Rate [Hz/m^{2}]");
    G1[0]->GetXaxis()->SetTitle("Voltage [V]");
    G1[0]->Draw("AP");
    
    //leg3->AddEntry(G1[0],"233 [K]","lp");
    
    //-110 C
    G1[4]->GetYaxis()->SetRangeUser(1e-3,1e9);
    G1[4]->GetYaxis()->SetLimits(1e-3,1e9);
    G1[4]->GetXaxis()->SetRangeUser(44,57);
    G1[4]->GetXaxis()->SetLimits(44,57);
    G1[4]->Draw("samep");
    
    //leg3->AddEntry(G1[4],"163 [K]","lp");
    
    //Enahnced contribution
    TF1 *fun_5=new TF1("fun_5",Fit_40_enhanced,45,57,21);
    
    fun_5->SetParameters(currentPar);
    fun_5->SetLineColor(Colors(0));
    fun_5->SetLineStyle(1);
    fun_5->SetLineWidth(2);
    fun_5->Draw("same,l");
    
    leg3->AddEntry(fun_5,"233 [K] SRH Enh.","lp");
    
    //Standard contribution
    TF1 *fun_6=new TF1("fun_5",Fit_40_standard,45,57,21);
    
    fun_6->SetParameters(currentPar);
    fun_6->SetLineColor(Colors(1));
    fun_6->SetLineStyle(1);
    fun_6->SetLineWidth(2);
    fun_6->Draw("same,l");
    
    leg3->AddEntry(fun_6,"233 [K] SRH","lp");

    //Sum
    TF1 *fun_7=new TF1("fun_0",Fit_40,45,57,21);
    
    fun_7->SetParameters(currentPar);
    fun_7->SetLineColor(Colors(0));
    fun_7->SetLineStyle(1);
    fun_7->SetLineWidth(2);
    
    //     mg1->Add(new TGraph(fun_1));
    //fun_7->Draw("same");
    
    //leg3->AddEntry(fun_7,"233 [K]","lp");
    
    
    //Enahnced contribution
    TF1 *fun_8=new TF1("fun_5",Fit_110_enhanced,44,57,21);
    
    fun_8->SetParameters(currentPar);
    fun_8->SetLineColor(Colors(4));
    fun_8->SetLineStyle(1);
    fun_8->SetLineWidth(2);
    fun_8->Draw("same,l");
    
    leg3->AddEntry(fun_8,"163 [K] SRH Enh.","lp");
    
    //Standard contribution
    TF1 *fun_9=new TF1("fun_5",Fit_110_standard,44,57,21);
    
    fun_9->SetParameters(currentPar);
    fun_9->SetLineColor(Colors(5));
    fun_9->SetLineStyle(1);
    fun_9->SetLineWidth(2);
    fun_9->Draw("same,l");

    leg3->AddEntry(fun_9,"163 [K] SRH","lp");
    
    //Sum
    TF1 *fun_10=new TF1("fun_0",Fit_110,44,57,21);
    
    fun_10->SetParameters(currentPar);
    fun_10->SetLineColor(Colors(4));
    fun_10->SetLineStyle(1);
    fun_10->SetLineWidth(2);
    
    //     mg1->Add(new TGraph(fun_1));
    //fun_10->Draw("same");
    
    //leg3->AddEntry(fun_10,"163 [K]","lp");
    
    leg3->Draw("same");
    

 
        
    std::cout<<" MINUIT FIT ......"<<std::endl;
    
    double chi2, edm, errdef;
    int nvpar, nparx;
    mini->GetStats(chi2,edm,errdef,nvpar,nparx);
    
    std::cout<<" chi2 "<<chi2<<" nvpar "<<nvpar<<" nparx "<<nparx<<std::endl;
    
    
    //Now compute total number of points you are fitting
    
    double Total_number_points=0;
    
    
    for(int ll=0;ll<NFILE;ll++) Total_number_points=Total_number_points+G1[ll]->GetN();
    
    
    std::cout<<" Total Number of Point is "<<Total_number_points<<std::endl;
    
    int ndf = Total_number_points-nvpar;
    
    std::cout<<" Residual Chi Square "<<chi2/ndf<<std::endl;
    
}



