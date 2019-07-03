#ifndef CONVECTED_DUBINS_PATH_H
#define CONVECTED_DUBINS_PATH_H

#include<vector>
#include<string>

#include<ConvectedDubins_Problem.h>

namespace ConvectedDubins {


// The RobustDubins::Path object can be used on its own (e.g., to draw a Dubins
// path with known parameters, or it can be the output of 
// the RobustDubins::Solver

enum pathClass { BSB, BBB, NO_CLASS };

enum pathType { LSL, RSR, LSR, RSL, LRL, RLR, NO_TYPE };

enum pathStatus { INFEASIBLE, FEASIBLE, OPTIMAL, NO_STATUS };

class Path {

private: 

  // Note: BSB paths are defined with
  // 1st segment: trochoid with t \in [0, tAA)
  // 2nd segment: straight line with t \in [tA, tBeta)
  // 3rd segment: trochoid with t\in [tBeta, T]
  pathClass pc_;
  pathType pt_;
  pathStatus ps_;   

  // to avoid storing redundant information, the ConvectedDubins::Path can 
  // reference an already defined problem 
  ConvectedDubins::Problem * prob_;
  ConvectedDubins::Problem probGen_; // used only when a problem is not supplied directly  
  bool probSuppliedFlag;

  // segment A: trochoid 
  double tA_; // time duration
  double dirA_; // direction 

  // segment B: straight, or trochoid 
  double tBeta_; // elapsed time to start trochoid B  

  // segment C: trochoid 
  double T_; // elapsed total time 
  double dirB_; // direction 
 
  // derived quantities: elapsed time from init. condition 
//  double phitA_; // phase
//  double phitB_; // phase

  // plotting
  double dt_; // time-step for defining state history
  vd x_, y_, h_, t_;  // state history in the inertial frame
  bool pathHistoryComputedFlag_;
//  vd xFi_, yFi_, hFi_;  // state history in inertial frame
 

  // constants 
  double t2pi_;
  void initialize();


  // private set functions
	void set_xInitial( double xInitial );
	void set_yInitial( double yInitial );
	void set_hInitialRad( double hInitialRad );

  // intermediate states at junction of extremals
	double xFinal_; // Final x position, default=0 (LU)
	double yFinal_; // Final y position, default=0 (LU)
	double hFinalRad_; // Final heading angle, default=0 (rad)
  double xBinit_, yBinit_, hBinitRad_;
  double xBfinal_, yBfinal_, hBfinalRad_;
  bool connectedStatesComputedFlag_;

  void pathXYH( const double & t, double & x, double & y, double & h);
  void printInputType( bool userInputFlag );
public: 
  // constructor and destruct
	Path();
  Path( ConvectedDubins::Problem * prob );
	virtual ~Path(){};

  void set_prob( ConvectedDubins::Problem * prob  );
  void set_windMagXY( double windMagX, double windMagY);
  void set_windMagDir( double windMagnitude, double windDirRad );
  void set_vehicleProperties( double minTurningRadius, double nominalSpeed );
  void set_stateInitial( const vd & stateInitial );
  void set_stateInitial( double xInitial, double yInitial, double hInitialRad );

  void set_pathParamsBSB( double dirA, double dirB, double tA, double tBeta, double T);
  void set_pathParamsBBB( double dirA, double tA, double tBeta, double T);
  void set_pathStatusOptimal();
  void set_pathStatusInfeasible();
  void set_pathType( ConvectedDubins::pathType pt){ pt_ = pt;};

  void print();  
  bool isDefined();
  ConvectedDubins::pathStatus get_pathStatus(){return ps_;};
  ConvectedDubins::pathClass get_pathClass(){return pc_;};
  ConvectedDubins::pathType get_pathTyped(){return pt_;};
  void computeConnectingStates();
  void computePathHistory();
  void printWaypoints();
  void plotXY(const std::string & fileName, const int & figureNumber);
  double get_T(){return T_;};


  void get_waypoints( vd & x, vd & y, vd & hRad );
  void get_waypoints_XYEastNorth( vd & x, vd & y, vd & hRad );

  double get_xFinal(){return xFinal_;};
  double get_yFinal(){return yFinal_;};
  double get_hFinalRad(){return hFinalRad_;};

}; // class

} // namespace
#endif
