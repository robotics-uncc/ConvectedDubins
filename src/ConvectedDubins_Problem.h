#ifndef CONVECTED_DUBINS_PROBLEM_H
#define CONVECTED_DUBINS_PROBLEM_H

#include<vector>

typedef std::vector<double> vd;

namespace ConvectedDubins {

// The ConvectedDubins::Problem object is mainly used to define a problem as an 
// input for the ConvectedDubins::Solver 

class Problem {

private:

  // Note: all angles are measured in radians, CW from x-axis which points north

	// optional input: final state 
	double xInitial_; // default=0
	double yInitial_; // default=0
	double hInitialRad_; // default=0
  bool startPtInputFlag_; // flags if user inputs a non-default init. cond.

  // required inputs 
	double xFinal_; // Final x position, default=0 (LU)
	double yFinal_; // Final y position, default=0 (LU)
	double hFinalRad_; // Final heading angle, default=0 (rad)
  bool endPtInputFlag_; // flag determining if problem is defined 

  // optional input: vehicle properties
	double R_; // minimum turning radius, default=1 (LU)
  double v_; // vehicle speed, default=1 (LU/s)
  double omega_; // derived quantity: vehicle maximum turn rate (rad/s)
  bool vehPropInputFlag_; // flags if user inputs a non-default radius

  // environment properties
  double w_; // wind magnitude, default=0.5(LU/s))
  double hwRad_; // wind direction, default=PI/2 (rad)
  double wx_; // wind x-component
  double wy_; // wind y-component
  bool envPropInputFlag_;

  // called by the constructor 
  void initialize();

  // private set functions
	void set_xInitial( double xInitial );
	void set_yInitial( double yInitial );
	void set_hInitialRad( double hInitialRad );

	void set_xFinal( double xFinal ); 
	void set_yFinal( double yFinal );
	void set_hFinalRad( double hFinalRad );

  void printInputType( bool userInputFlag );

public:

	// constructor and destructor
  Problem( );
	virtual ~Problem(){};

	// set functions
	void set_stateInitial( double xInitial, 
                         double yInitial, 
                         double hInitialRad );

	void set_stateInitial_XYEastNorth( double xInitial, 
                                     double yInitial, 
                                     double hInitialRad );
	void set_stateInitial_XYEastNorthCourse( double xInitial, 
                                           double yInitial, 
                                           double chiInitialRad );
  
	void set_stateInitial( const vd & stateInitial );

	void set_stateFinal( double xFinal, 
                       double yFinal, 
                       double hFinalRad );
	void set_stateFinal_XYEastNorth( double xFinal, 
                                   double yFinal, 
                                   double hFinalRad );
	void set_stateFinal_XYEastNorthCourse( double xFinal, 
                                         double yFinal, 
                                         double chiFinalRad );

	void set_stateFinal( const vd & stateFinal );

  void set_vehicleProperties( double minTurningRadius, 
                              double nominalSpeed );

  void set_windMagDir( double windMagnitude, double windDirRad);
  void set_windMagXY( double windMagX,  double windMagY);

	// main functions
	void print(); // prints problem summary to screen 

	// get functions
  bool isDefined();
  bool isDefinedNoEndpoint();

	double xInitial(){return xInitial_;}; 
	double yInitial(){return yInitial_;};
	double hInitialRad(){return hInitialRad_;};
	vd stateInitial();

	double xFinal();
	double yFinal();
	double hFinalRad();
	vd stateFinal();

  double R(){return R_;};
  double v(){return v_;};
  double omega(){return omega_;};

  double w(){return w_;};
  double hwRad(){return hwRad_;};
  double wx(){return wx_;};
  double wy(){return wy_;};


}; // class

} // namespace

#endif


