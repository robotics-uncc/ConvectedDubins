#include<string>
#include<vector>
#include<cmath>
#include<iostream>

// custom
#include<ConvectedDubins_Problem.h>
#include<ConvectedDubins_Path.h>

int main() {

  ConvectedDubins::Path path;
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
  path.set_pathParamsBBB( dirA, tA, tBeta, T);



  double dirB;
  if ( defaultFlag == 1){
    path.set_stateInitial(x0,y0,h0Deg*M_PI/180.0);
    path.set_vehicleProperties(R, v);
    path.set_windMagDir(w, dirDeg*M_PI/180.0);
  }
  else if ( defaultFlag == 0){
    printf("Initial X-Pos : ");
    std::cin  >> x0; 
    printf("Initial Y-Pos : ");
    std::cin  >> y0; 
    printf("Initial Heading Deg : ");
    std::cin  >> h0Deg; 
    path.set_stateInitial(x0,y0,h0Deg*M_PI/180.0);

    printf("Turn Radius : ");
    std::cin  >> R;
    printf("Nom. Vehicle Speed : ");
    std::cin  >> v;
    path.set_vehicleProperties(R, v);

    printf("Wind Magnitude : ");
    std::cin  >> w;
    printf("Wind Direction Deg : ");
    std::cin  >> dirDeg;
    path.set_windMagDir(w, dirDeg*M_PI/180.0);

    int pathType; 
    printf(" Path Class? (0 - BSB, 1 - BBB): ");
    std::cin >> pathType;
    printf(" Initial Trochoid Dir (1 : CCW, -1 : CW) : ");
    std::cin >> dirA;
    printf(" First-to-Second Segment Switch Time : ");
    std::cin >> tA;
    if ( pathType == 0 ){
      double straight_time;
      printf(" Second-to-Third Segment Switch Time : ");
      std::cin >> tBeta;
      printf(" Final Trochoid Dir (1 : CCW, -1 : CW) : ");
      std::cin >> dirB;
      printf(" Terminal Time : ");
      std::cin >> T;
      path.set_pathParamsBSB( dirA, dirB, tA, tBeta, T);
    }
    else if (pathType == 1) {
      printf(" Second-to-Third Segment Switch Time : ");
      std::cin >> tBeta;
      printf(" Final Trochoid Dir (1 : CCW, -1 : CW) : ");
      std::cin >> dirB;
      printf(" Terminal Time : ");
      std::cin >> T;
      path.set_pathParamsBBB( dirA, tA, tBeta, T);
    }
    else { 
      printf("invalid entry.\n");
      return 0;
    }
  }
  else {
    printf("invalid entry.\n");
    return 0;
  }



  printf("plot_LSL_RSR(1.0)\n");
  path.print();
  path.computePathHistory();
  printf("printWaypoints()\n");
  //path.printWaypoints();
  path.plotXY("test.m",1);
  
 
	return 0;
}




