#include<string>
#include<vector>
#include<cmath>
#include<iostream>

// custom
#include<ConvectedDubins_Problem.h>
#include<ConvectedDubins_Path.h>
#include<ConvectedDubins_Solver.h>

int main() {

  printf("--------------------------------------------\n");
  printf("ConvectedDubins: Path Plotter\n");
  printf("--------------------------------------------\n");

  int defaultFlag;
  printf("Use default path? (0-no, 1-yes)");
  std::cin >> defaultFlag;


  // example Fig. 12, JGCD paper (Techy et al))
  double x0 = 0.0;
  double y0 = -200.0;
  double h0Deg = 0.0;

  double x1 = 100.0;
  double y1 = 100.0;
  double h1Deg = 90.0;

  double v = 20;
  double omega = 0.2832;
  double R = v/omega;
  double w = 5.0;
  double dirDeg = 0;
  double dirA = -1;
  double tA = 0.6869;
  double T = 22.4228;  
  double tBeta = 17.4451;
  double t2pi = 2.0*M_PI/(v/R);

  double straight_time, trochoidC_dir;

  ConvectedDubins::Problem prob;
  if ( defaultFlag == 1){
    prob.set_stateInitial(x0,y0,h0Deg*M_PI/180.0);
    prob.set_stateFinal(x1,y1,h1Deg*M_PI/180.0);
    prob.set_vehicleProperties(R, v);
    prob.set_windMagDir(w, dirDeg*M_PI/180.0);
  }
  else if ( defaultFlag == 0){
    printf("Initial X-Pos : ");
    std::cin  >> x0; 
    printf("Initial Y-Pos : ");
    std::cin  >> y0; 
    printf("Initial Heading Deg : ");
    std::cin  >> h0Deg; 
    prob.set_stateInitial(x0,y0,h0Deg*M_PI/180.0);

    printf("Final X-Pos : ");
    std::cin  >> x1; 
    printf("Final Y-Pos : ");
    std::cin  >> y1; 
    printf("Final Heading Deg : ");
    std::cin  >> h1Deg; 
    prob.set_stateFinal(x1, y1, h1Deg*M_PI/180.0);

    printf("Turn Radius : ");
    std::cin  >> R;
    printf("Nom. Vehicle Speed : ");
    std::cin  >> v;
    prob.set_vehicleProperties(R, v);

    printf("Wind Magnitude : ");
    std::cin  >> w;
    printf("Wind Direction Deg : ");
    std::cin  >> dirDeg;
    prob.set_windMagDir(w, dirDeg*M_PI/180.0);
  }
  else {
    printf("invalid entry.\n");
    return 0;
  }

  prob.print();
  printf(" set prob. \n");
  ConvectedDubins::Path path( &prob );
  printf(" set params. \n");
  path.set_pathParamsBBB( dirA, tA, tBeta, T);
  path.print();
  path.computePathHistory();
  printf("printWaypoints()\n");
//  //path.printWaypoints();
  path.plotXY("test.m",1);
  
 
	return 0;
}




