#include<limits> 
#include<ConvectedDubins_Path.h>
#include<ConvectedDubins_Utils.h>
#include<stdexcept>
#include<MathTools.h>

ConvectedDubins::Path::Path(){
  initialize();
}

ConvectedDubins::Path::Path( ConvectedDubins::Problem * prob ){
  set_prob( prob );
  initialize();
}

void ConvectedDubins::Path::set_prob( ConvectedDubins::Problem * prob  ){
  #ifndef NDEBUG
  if ( !prob->isDefinedNoEndpoint() ){
    throw std::runtime_error("ConvectedDubins::Path::Path: Problem not defined sufficiently.");
  }
  #endif 
  prob_ = prob;
  probSuppliedFlag = true;  
}


void ConvectedDubins::Path::initialize(){
  tA_ = 0.0;
  tBeta_ = 0.0;
  T_ = 0.0;
  dt_ = 0.1;
  dirA_ = 0.0;
  dirB_ = 0.0;
  ps_ = ConvectedDubins::pathStatus::NO_STATUS;
  pt_ = ConvectedDubins::pathType::NO_TYPE;
  pc_ = ConvectedDubins::pathClass::NO_CLASS;
  pathHistoryComputedFlag_ = false;
  connectedStatesComputedFlag_ = false;
  if ( probSuppliedFlag == false ){
    prob_ = &probGen_;
  }
  
}

void ConvectedDubins::Path::set_pathParamsBSB( double dirA, double dirB,
                                               double tA, double tBeta, double T){
  dirA_ = dirA;
  dirB_ = dirB;
  tA_ = tA;
  tBeta_ = tBeta;
  T_ = T;
  pc_ = ConvectedDubins::pathClass::BSB;
  if ( dirA_ == -1 && dirB_ == -1){
    pt_ = ConvectedDubins::pathType::LSL;
  }
  else if ( dirA_ == 1 && dirB_ == 1){
    pt_ = ConvectedDubins::pathType::RSR;
  }
  else if ( dirA_ == -1 && dirB_ == 1){
    pt_ = ConvectedDubins::pathType::LSR;
  }
  else if ( dirA_ == 1 && dirB_ == -1){
    pt_ = ConvectedDubins::pathType::RSL;
  }
  else {
    throw std::runtime_error("ConvectedDubins::Path::set_pathParamsBSB : dirA and dirB invalid.");
  }
  ps_ = ConvectedDubins::pathStatus::FEASIBLE;
}

void ConvectedDubins::Path::set_pathParamsBBB( double dirA, double tA, double tBeta, double T){
  dirA_ = dirA;
  dirB_ = dirA;
  tA_ = tA;
  tBeta_ = tBeta;
  T_ = T;
  pc_ = ConvectedDubins::pathClass::BBB;
  if ( dirA_ == -1 ){
    pt_ = ConvectedDubins::pathType::LRL;
  }
  else if ( dirA_ == 1 ){
    pt_ = ConvectedDubins::pathType::RLR;
  }
  else {
    throw std::runtime_error("ConvectedDubins::Path::set_pathParamsBBB : dirA invalid.");
  }
  ps_ = ConvectedDubins::pathStatus::FEASIBLE;
}

void ConvectedDubins::Path::set_pathStatusOptimal(){
  ps_ = ConvectedDubins::pathStatus::OPTIMAL;
}

void ConvectedDubins::Path::set_pathStatusInfeasible(){
  ps_ = ConvectedDubins::pathStatus::INFEASIBLE;
}

void ConvectedDubins::Path::set_stateInitial( double xInitial, 
                                              double yInitial, 
                                              double hInitialRad ){
  prob_->set_stateInitial( xInitial, yInitial, hInitialRad );
}

void ConvectedDubins::Path::set_stateInitial( const vd & stateInitial ){
	prob_->set_stateInitial( stateInitial );
}


void ConvectedDubins::Path::set_vehicleProperties( double minTurningRadius, 
                                                   double nominalSpeed ){
  prob_->set_vehicleProperties( minTurningRadius, nominalSpeed );
}

void ConvectedDubins::Path::set_windMagDir( double windMagnitude, 
                                            double windDirRad ){
  prob_->set_windMagDir( windMagnitude, windDirRad );
}

void ConvectedDubins::Path::set_windMagXY( double windMagX,  
                                           double windMagY){
  prob_->set_windMagXY( windMagX, windMagY ); 
}


void ConvectedDubins::Path::computeConnectingStates(){

  // state after first trochoid straight
  xBinit_ = prob_->xInitial() + trochoid_delx( prob_->v(), prob_->omega(), prob_->wx(), dirA_, 
                                               tA_, prob_->hInitialRad() );
  yBinit_ = prob_->yInitial() + trochoid_dely( prob_->v(), prob_->omega(), prob_->wy(), dirA_, 
                                               tA_, prob_->hInitialRad() );
  hBinitRad_ = fmod( prob_->hInitialRad() + trochoid_delh( prob_->omega(), dirA_, tA_ ) , 
                     2.0*M_PI );
  if ( hBinitRad_ < 0 ){
    hBinitRad_ = hBinitRad_ + 2.0*M_PI;
  }
  switch( pc_ ){
    // state after second segment (straight)
    case ConvectedDubins::pathClass::BSB: {
      xBfinal_ = xBinit_ + straight_delx( prob_->v(), prob_->wx(), tBeta_ - tA_ , hBinitRad_ );
      yBfinal_ = yBinit_ + straight_dely( prob_->v(), prob_->wy(), tBeta_ - tA_ , hBinitRad_ );
      hBfinalRad_ = hBinitRad_;
      break;
    }
    // state after second segment (trochoid)
    case ConvectedDubins::pathClass::BBB: {
      xBfinal_ = xBinit_ + trochoid_delx( prob_->v(), prob_->omega(), prob_->wx(), -dirA_, 
                                          tBeta_ - tA_ , hBinitRad_ );
      yBfinal_ = yBinit_ + trochoid_dely( prob_->v(), prob_->omega(), prob_->wy(), -dirA_, 
                                          tBeta_ - tA_ , hBinitRad_ );
      hBfinalRad_ = fmod( hBinitRad_ + trochoid_delh( prob_->omega(), -dirA_, tBeta_ - tA_ ) , 
                          2.0*M_PI ); 
      if ( hBfinalRad_ < 0 ){
        hBfinalRad_ = hBfinalRad_ + 2.0*M_PI;
      }
      break;
    }
  }

  // final state 
  xFinal_ = xBfinal_ + trochoid_delx( prob_->v(), prob_->omega(), prob_->wx(), dirB_, 
                                      T_ - tBeta_ , hBfinalRad_ );
  yFinal_ = yBfinal_ + trochoid_dely( prob_->v(), prob_->omega(), prob_->wy(), dirB_, 
                                      T_ - tBeta_ , hBfinalRad_ );
  hFinalRad_ = fmod( hBfinalRad_ + trochoid_delh( prob_->omega(), dirB_, T_ - tBeta_ ) , 2.0*M_PI ); 
  if ( hFinalRad_ < 0 ){
    hFinalRad_ = hFinalRad_ + 2.0*M_PI;
  }  
//  printf(" path class : %i \n", pc_);
//  printf(" (xBinit, yBinit, hBinit) = (%3.3f, %3.3f, %3.3f) \n", xBinit_, yBinit_, hBinitRad_);
//  printf(" (xBfinal, yBfinal, hBfinal) = (%3.3f, %3.3f, %3.3f) \n", xBfinal_, yBfinal_, hBfinalRad_);
//  printf(" (xFinal, yFinal, hFinal) = (%3.3f, %3.3f, %3.3f) \n", xFinal_, yFinal_, hFinalRad_);

  connectedStatesComputedFlag_ = true;
}

bool ConvectedDubins::Path::isDefined(){
  if ( ps_ == ConvectedDubins::pathStatus::FEASIBLE ||
       ps_ == ConvectedDubins::pathStatus::OPTIMAL ){
    return true;
  }
  else {
    return false;
  }
}

void ConvectedDubins::Path::pathXYH( const double & t, double & x, double & y, double & h){
  if ( isDefined() ){
    if ( !connectedStatesComputedFlag_ ){
      computeConnectingStates();
    }
    if ( t <= tA_ && t >= 0 ){ // first segment: trochoid 
      x = prob_->xInitial() + trochoid_delx( prob_->v(), prob_->omega(), prob_->wx(), dirA_, 
                                             t, prob_->hInitialRad() );
      y = prob_->yInitial() + trochoid_dely( prob_->v(), prob_->omega(), prob_->wy(), dirA_, 
                                             t, prob_->hInitialRad() );
      h = prob_->hInitialRad() + trochoid_delh( prob_->omega(), dirA_, t );
      //printf("Seg. A: (t,x,y,h) = (%3.3f,%3.3f,%3.3f,%3.3f) \n",t,x,y,h);
    }
    else if ( t <= tBeta_ ){
      switch ( pc_ ){
        case ConvectedDubins::pathClass::BSB: { // second segment: straight line 
          x = xBinit_ + straight_delx( prob_->v(), prob_->wx(), t - tA_ , hBinitRad_ );
          y = yBinit_ + straight_dely( prob_->v(), prob_->wy(), t - tA_ , hBinitRad_ );
          h = hBinitRad_;
          //printf("Seg. B / BSB: (t,x,y,h) = (%3.3f,%3.3f,%3.3f,%3.3f) \n",t,x,y,h);
          break;
        }
        case ConvectedDubins::pathClass::BBB: { // second segment: trochoid 
          x = xBinit_ + trochoid_delx( prob_->v(), prob_->omega(), prob_->wx(), -dirA_, 
                                       t - tA_ , hBinitRad_ );
          y = yBinit_ + trochoid_dely( prob_->v(), prob_->omega(), prob_->wy(), -dirA_, 
                                       t - tA_ , hBinitRad_ );
          h = hBinitRad_ + trochoid_delh( prob_->omega(), -dirA_, t - tA_ );
          //printf("Seg. B / BBB: (t,x,y,h) = (%3.3f,%3.3f,%3.3f,%3.3f) \n",t,x,y,h);
          break;
        }  
        
      }
    }
    else if ( t <= T_){ // third segment: trochoid 
        x = xBfinal_ + trochoid_delx( prob_->v(), prob_->omega(), prob_->wx(), dirB_, 
                                      t - tBeta_ , hBfinalRad_ );
        y = yBfinal_ + trochoid_dely( prob_->v(), prob_->omega(), prob_->wy(), dirB_, 
                                      t - tBeta_ , hBfinalRad_ );
        h = hBfinalRad_ + trochoid_delh( prob_->omega(), dirB_, t - tBeta_ );
        //printf("Seg. C: (t,x,y,h) = (%3.3f,%3.3f,%3.3f,%3.3f) \n",t,x,y,h);
    }
    else {
       throw std::runtime_error("ConvectedDubins::Path::pathXYH: incorrect time.");
    }
  }
  else {
    throw std::runtime_error("ConvectedDubins::Path::pathXYH: path not defined.");
  }
  h = fmod( h, 2.0*M_PI );
  if ( h < 0 ){
    h = h + 2.0*M_PI;
  }
}

void ConvectedDubins::Path::computePathHistory(){
  if ( isDefined() ){
    int N = (int)floor(T_/dt_);
    //printf(" N = %i, T = %3.3f, dt = %3.3f \n",N,T_, dt_);
    x_.resize(N);
    y_.resize(N);
    h_.resize(N);
    t_.resize(N);
    double x,y,h;
    for (int i=0; i < N-1; i++){
      t_[i] = i*dt_;
      pathXYH( t_[i], x, y, h);
      x_[i] = x;
      y_[i] = y;
      h_[i] = h;
    }
    // force endpoint 
    t_[N-1] = T_;
    pathXYH( t_[N-1], x, y, h);
    x_[N-1] = x;
    y_[N-1] = y;
    h_[N-1] = h;
    pathHistoryComputedFlag_ = true;
  }
  else {
    std::runtime_error("ConvectedDubins::Path::computePathHistory(): path is undefined.");
  }
}

void ConvectedDubins::Path::get_waypoints( vd & x, vd & y, vd & hRad ){
  if ( !pathHistoryComputedFlag_ ){
    computePathHistory();
  }
  if ( isDefined() ){
    x = x_;
    y = y_;
    hRad = h_;
  }
}

void ConvectedDubins::Path::get_waypoints_XYEastNorth( vd & x, vd & y, vd & hRad ){
  if ( !pathHistoryComputedFlag_ ){
    computePathHistory();
  }
  if ( isDefined() ){
    x = y_;
    y = x_;
    vd h_XYEastNorth;
    for ( int i = 0; i < h_.size(); i++){
      h_XYEastNorth.push_back( fmod( -h_[i] + M_PI/2.0 , 2.0*M_PI ) );
    }
    hRad = h_XYEastNorth;
  }
}


void ConvectedDubins::Path::plotXY(const std::string & fileName, 
                                   const int & figureNumber){
  if ( (ps_ == FEASIBLE || ps_ == OPTIMAL) && pathHistoryComputedFlag_ ){
    // write path plotting commands
	  MathTools::writeVectorData(fileName, x_, "x");
	  MathTools::writeVectorData(fileName, y_, "y");
    std::string pathStyle = "k--";
    if ( ps_ == ConvectedDubins::pathStatus::OPTIMAL ){
      pathStyle = "b-";
    }
	  MathTools::plot1DCurve(fileName, figureNumber, "y", "x", pathStyle);
    MathTools::axisProperty(fileName, "equal");
  }
  else {
    printf("ConvectedDubins::Path::plotXY: cannot plot because path is undefined.\n");
  }
}

void ConvectedDubins::Path::printWaypoints(){
  if ( pathHistoryComputedFlag_ ){
    for (int i=0; i < x_.size(); i++){
      printf(" (x%i,y%i) = (%3.3f, %3.3f) \n",i,i,x_[i],y_[i]);
    }
  }
  else {
    printf("ConvectedDubins::Path::printWaypoints(): cannot print because path is undefined.");
  }
}

void ConvectedDubins::Path::print(){
  printf("Convected Dubins Path \n");
    switch ( pt_ ){
      case ConvectedDubins::pathType::NO_TYPE: printf("\tType: undefined \n"); break;
      case ConvectedDubins::pathType::LSL: printf("\tType: LSL \n"); break;
      case ConvectedDubins::pathType::RSR: printf("\tType: RSR \n"); break;
      case ConvectedDubins::pathType::LSR: printf("\tType: LSR \n"); break;
      case ConvectedDubins::pathType::RSL: printf("\tType: RSL \n"); break;
      case ConvectedDubins::pathType::LRL: printf("\tType: LRL \n"); break;
      case ConvectedDubins::pathType::RLR: printf("\tType: RLR \n"); break;
    }
    switch ( ps_ ){
      case pathStatus::NO_STATUS: printf("\tSolution Status: undefied \n"); break;
      case pathStatus::INFEASIBLE: printf("\tSolution Status: infeasible \n"); break;
      case pathStatus::FEASIBLE: printf("\tSolution Status: feasible \n"); break;
      case pathStatus::OPTIMAL: printf("\tSolution Status: optimal \n"); break;
    }
    if ( ps_ == FEASIBLE || ps_ == OPTIMAL ){
      printf("\tSeg. A Endpoint Elasped Time (tA) : %3.3f \n",tA_);
      printf("\tSeg. B Endpoint Elasped Time (tBeta) : %3.3f \n",tBeta_);
      printf("\tSeg. C Endpoint/Total Elasped Time (T) : %3.3f \n",T_);
    }
    if ( (ps_ == FEASIBLE || ps_ == OPTIMAL) && connectedStatesComputedFlag_ ){
      printf("\tInitial State: \t(x,y,hdeg) = (%3.3f,\t %3.3f,\t %3.3f) \n",prob_->xInitial(), 
                                                      prob_->yInitial(), 
                                                      prob_->hInitialRad()*180.0/M_PI );
      printf("\tSeg. A Endpoint : \t(x,y,hdeg) = (%3.3f,\t %3.3f,\t %3.3f) \n", xBinit_, 
                                                                         yBinit_, 
                                                                   hBfinalRad_*180.0/M_PI);
      printf("\tSeg. B Endpoint : \t(x,y,hdeg) = (%3.3f,\t %3.3f,\t %3.3f) \n", xBfinal_, 
                                                                         yBfinal_, 
                                                                    hBinitRad_*180.0/M_PI);
      printf("\tFinal State: \t\t(x,y,hdeg) = (%3.3f,\t %3.3f,\t %3.3f) \n", xFinal_, yFinal_, 
                                                                hFinalRad_*180.0/M_PI);
    }
  
  
}

