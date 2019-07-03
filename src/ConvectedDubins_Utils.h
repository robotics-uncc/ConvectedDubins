#ifndef CONVECTED_DUBINS_UTILS_H
#define CONVECTED_DUBINS_UTILS_H

#include<ConvectedDubins_Path.h>
#include<vector>
#include<string>
#include<assert.h>

typedef std::vector<double> vd;
typedef std::vector<int> vi;


namespace ConvectedDubins {


// v - flow-relative speed
// w - disturbance (current/wind) speed
// hwRad - direction of disturbance, measured CCW from x axis poiniting East
// chiRad - course angle of vehicle's instantaenous motion 
double headingAngleFromInertialCourse( double v, double w, double hWRad, 
                                       double chiRad);

double speedOverGround( const double & va, 
                        const double & vw, 
                        const double & hRad);

double courseAngleRad( const double & va, 
                       const double & vw, 
                       const double & hRad);

//void convertXYTrochoidToInertial( const vd & xt, const vd & yt, 
//                                  const double & hwRad, 
//                                  vd & xN, vd & yE);

void convertInertialtoTrochoid( double xn, double ye,
                                double hwRad,
                                double & xt, double & yt );


ConvectedDubins::pathClass getClassFromType( ConvectedDubins::pathType pt );

double trochoid_delx ( double v, double omega, double wx, 
                       double dir, double delt, double hInitRad );
double trochoid_delx ( ConvectedDubins::Problem * prob, 
                       double dir, double delt, double hInitRad );


double trochoid_dely ( double v, double omega, double wy, 
                       double dir, double delt, double hInitRad );
double trochoid_dely ( ConvectedDubins::Problem * prob, 
                       double dir, double delt, double hInitRad );

double trochoid_delh ( double omega, 
                       double dir, double delt );
double trochoid_delh ( ConvectedDubins::Problem * prob, 
                       double dir, double delt );

double straight_delx ( double v, double wx, 
                       double delt, double hInitRad );
double straight_delx ( ConvectedDubins::Problem * prob, 
                       double delt, double hInitRad );

double straight_dely ( double v, double wy, 
                       double delt, double hInitRad );
double straight_dely ( ConvectedDubins::Problem * prob, 
                       double delt, double hInitRad );


// tA is the end of first segment 
// tBeta is the end of second segment
// T is the end of third segment (total elapsed time))
void computeEndpoint( ConvectedDubins::Problem * prob, 
                      ConvectedDubins::pathType pt,
                      double dirA, double dirC,
                      double tA, double tBeta, double T,
                      double & xFinal, double & yFinal, 
                      double & hFinalRad);

// like fmod, but ensures returned value is positive between [0, mod]]
double fmodPos( double val, double mod );


} // namespace

#endif
