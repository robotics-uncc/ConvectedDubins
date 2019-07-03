#include<string>
#include<vector>
#include<cmath>
#include<iostream>

// custom
#include<ConvectedDubins_Problem.h>
#include<ConvectedDubins_Path.h>
#include<ConvectedDubins_Solver.h>
#include<MathTools.h>
int main() {

  ConvectedDubins::Path path;
  printf("--------------------------------------------\n");
  printf("ConvectedDubins: Path Plotter\n");
  printf("--------------------------------------------\n");

  int defaultFlag;
  printf("Use default path? (0-no, 1-yes)");
  std::cin >> defaultFlag;


  double x0 = 0.0;
  double y0 = -200.0;
  double h0Deg = 0.0;

  double x1 = 0.0;
  double y1 = -180.0;
  double h1Deg = 180; 

  double v = 20;
  double R = v/0.2832;
  double w = 5.0;
  double dirDeg = 0;

//  double trochoidA_time = 0.6869;
//  double trochoidA_dir = -1;
//  double trochoidB_time = 17.4451;
//  double trochoidC_time = 22.4228 - trochoidA_time - trochoidB_time;


  double t2pi = 2.0*M_PI/(v/R);
  double straight_time, trochoidC_dir;

  ConvectedDubins::Problem prob;

//  path.set_pathParamsBBB( trochoidA_time, trochoidA_dir, trochoidB_time, 
//                          trochoidC_time);


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
  printf("solving... \n");
  ConvectedDubins::Solver solver( &prob );
  std::string fileName = "test.m";
  solver.computeConnectingStates();
  solver.printSolns();
  solver.plotSolns(fileName);
  MathTools::runOctaveScript(fileName);

  
 
	return 0;
}




