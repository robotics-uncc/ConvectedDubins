#include<ConvectedDubins_Utils.h>
#include<ConvectedDubins_Problem.h>
#include<cmath>
#include<stdexcept>


double ConvectedDubins::headingAngleFromInertialCourse( double v, double w, double hWRad, 
                                                        double chiRad){
  return fmodPos( chiRad - asin( w/v*sin(hWRad-chiRad) ), 2.0*M_PI );
}

double ConvectedDubins::speedOverGround( const double & va, 
                                         const double & vw, 
                                         const double & hRad){
  double term = va*va + vw*vw + 2*va*vw*cos(hRad);
  return sqrt(term*term);
}

// In the trochoid frame!
double ConvectedDubins::courseAngleRad( const double & va, 
                                        const double & vw, 
                                        const double & hRad){
  double num = va*sin(hRad);
  double denom = va*cos(hRad) + vw;
  return fmod( atan2(num, denom) , 2.0*M_PI);
}


double ConvectedDubins::fmodPos( double val, double mod ){
  while ( val < 0 ){
    val = val + mod;
  }
  return fmod(val,mod);
}

void ConvectedDubins::convertInertialtoTrochoid( double xn, double ye,
                                                 double hwRad, 
                                                 double & xt, double & yt ){
    //printf(" (xn, ye, hwRad) = (%3.3f, %3.3f, %3.3f) \n",xn,ye,hwRad);
    xt =  cos(hwRad)*xn + sin(hwRad)*ye;
    yt = -sin(hwRad)*xn + cos(hwRad)*ye; 
    //printf(" (xt,yt) = (%3.3f, %3.3f) \n", xt, yt);
}


double ConvectedDubins::trochoid_delx( double v, double omega, double wx, 
                                       double dir, double delt, double hInitRad ){
  return v/(dir*omega)*( sin(dir*omega*delt+hInitRad) - sin(hInitRad) )+wx*delt;
}

double ConvectedDubins::trochoid_delx( ConvectedDubins::Problem * prob, 
                                       double dir, double delt, double hInitRad ){
  return prob->v()/(dir*prob->omega())*( sin(dir*prob->omega()*delt+hInitRad) 
                                         - sin(hInitRad) ) + prob->wx()*delt;
}

double ConvectedDubins::trochoid_dely( double v, double omega, double wy,
                                       double dir, double delt, double hInitRad ){
  return v/(dir*omega)*(-cos(dir*omega*delt+hInitRad) + cos(hInitRad) )+wy*delt;
}

double ConvectedDubins::trochoid_dely( ConvectedDubins::Problem * prob, 
                                       double dir, double delt, double hInitRad ){
  return prob->v()/(dir*prob->omega())*( -cos(dir*prob->omega()*delt+hInitRad) + cos(hInitRad) )
                  + prob->wy()*delt;
}

double ConvectedDubins::trochoid_delh( double omega, 
                                       double dir, double delt ){
  return dir*(omega)*delt;
}

double ConvectedDubins::trochoid_delh( ConvectedDubins::Problem * prob, 
                                       double dir, double delt ){
  return dir*( prob->omega() )*delt;
}

double ConvectedDubins::straight_delx( double v, double wx, 
                                       double delt, double hInitRad ){
  return ( v*cos(hInitRad) + wx )*delt;
}

double ConvectedDubins::straight_delx( ConvectedDubins::Problem * prob, 
                                       double delt, double hInitRad ){
  return ( prob->v()*cos(hInitRad) + prob->wx() )*delt;
}

double ConvectedDubins::straight_dely( double v, double wy, 
                                       double delt, double hInitRad ){
  return ( v*sin(hInitRad) + wy )*delt;
}

double ConvectedDubins::straight_dely( ConvectedDubins::Problem * prob, 
                                       double delt, double hInitRad ){
  return ( prob->v()*sin(hInitRad) + prob->wy() )*delt;
}

ConvectedDubins::pathClass ConvectedDubins::getClassFromType( 
                                                  ConvectedDubins::pathType pt){
  switch (pt){
    case ConvectedDubins::pathType::LSL: return ConvectedDubins::pathClass::BSB;
    case ConvectedDubins::pathType::LSR: return ConvectedDubins::pathClass::BSB;
    case ConvectedDubins::pathType::RSL: return ConvectedDubins::pathClass::BSB; 
    case ConvectedDubins::pathType::RSR: return ConvectedDubins::pathClass::BSB;
    case ConvectedDubins::pathType::LRL: return ConvectedDubins::pathClass::BBB;
    case ConvectedDubins::pathType::RLR: return ConvectedDubins::pathClass::BBB; 
  }

}

void ConvectedDubins::computeEndpoint( ConvectedDubins::Problem * prob, 
                                       ConvectedDubins::pathType pt,
                                       double dirA, double dirC,
                                       double tA, double tBeta, double T,
                                       double & xFinal, double & yFinal, 
                                       double & hFinalRad){

  // state after first trochoid straight
  double xBinit = prob->xInitial() + trochoid_delx( prob, dirA, tA, prob->hInitialRad() );
  double yBinit = prob->yInitial() + trochoid_dely( prob, dirA, tA, prob->hInitialRad() );
  double hBinit = prob->hInitialRad() + trochoid_delh( prob->omega(), dirA, tA ); 
  ConvectedDubins::pathClass pc = ConvectedDubins::getClassFromType(pt);


  double xBfinal, yBfinal, hBfinal;
  switch( pc ){
    // state after second segment (straight)
    case ConvectedDubins::pathClass::BSB: {
      xBfinal = xBinit + straight_delx( prob, tBeta - tA, hBinit );
      yBfinal = yBinit + straight_dely( prob, tBeta - tA, hBinit );
      hBfinal = hBinit;
      break;
    }
    // state after second segment (trochoid)
    case ConvectedDubins::pathClass::BBB: {
      xBfinal = xBinit + trochoid_delx( prob, -dirA, tBeta - tA , hBinit);
      yBfinal = yBinit + trochoid_dely( prob, -dirA, tBeta - tA , hBinit);
      hBfinal = hBinit + trochoid_delh( prob, -dirA, tBeta - tA ); 
      break;
    }
  }

  // final state 
  xFinal = xBfinal + trochoid_delx( prob, dirC, T-tBeta , hBfinal );
  yFinal = yBfinal + trochoid_dely( prob, dirC, T-tBeta , hBfinal );
  hFinalRad = hBfinal + trochoid_delh( prob, dirC, T-tBeta ); 
  hFinalRad = fmod(hFinalRad, 2.0*M_PI);
  if ( hFinalRad < 0 ){
    hFinalRad = hFinalRad + 2.0*M_PI;
  }
  
  //printf(" (xFinal, yFinal, hFinal) = (%3.3f, %3.3f, %3.3f) \n", xFinal, yFinal, hFinalRad);
}

