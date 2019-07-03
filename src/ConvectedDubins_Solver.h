#ifndef CONVECTED_DUBINS_SOLVER_H
#define CONVECTED_DUBINS_SOLVER_H

#include<ConvectedDubins_Problem.h>
#include<ConvectedDubins_Path.h>

namespace ConvectedDubins {

class Solver {

private: 
  ConvectedDubins::Problem * prob_;

  // initial and final state in the trochoidal frame 
  double x0_, y0_;
  double x1_, y1_;

  // intermediate variables used by solver
  double t2pi_; 

  // paths 
  ConvectedDubins::Path LSL_;
  ConvectedDubins::Path RSR_;
  ConvectedDubins::Path LSR_;
  ConvectedDubins::Path RSL_;
  ConvectedDubins::Path RLR_;
  ConvectedDubins::Path LRL_;
  
  ConvectedDubins::Path * optPath_;

  bool checkConditionsBSB( double d1, double d2, double tA, double tB, 
                           double xt10, double yt10, double phit1,
                           double xt20, double yt20, double phit2,
                           ConvectedDubins::pathType pt  );
  bool checkConditionsBBB( double d1, double tA, double tB, double T,
                           double xt10, double yt10, double phit1,
                           ConvectedDubins::pathType pt );

  void solveBSB( ConvectedDubins::pathType pt );
  void solveBBB( ConvectedDubins::pathType pt );
  double fixedPointBSB( double p0, double d1, double d2, double k,
                        double xt10, double yt10, double phit1,
                        double xt20, double yt20, double phit2 );
  bool fixedPointBBB( std::vector<double> pvecInit, double d1, double k, 
                      std::vector<double> & pvecOut,
                      double xt10, double yt10, double phit1 );
  void setOptimal();

public: 
  // constructor and destruct
	Solver( ConvectedDubins::Problem * prob );
	virtual ~Solver(){};

  ConvectedDubins::Path get_path( ConvectedDubins::pathType type );
  ConvectedDubins::Path * get_optPath(){return optPath_;};
  void get_optimalWaypoints( vd & x, vd & y, vd & hRad );
  void get_optimalWaypoints_XYEastNorth( vd & x, vd & y, vd & hRad );
  double get_optCost(){ return optPath_->get_T(); };

  void computeConnectingStates();
  void plotSolns(std::string fileName);
  void printSolns();



}; // class




} // namespace


#endif
