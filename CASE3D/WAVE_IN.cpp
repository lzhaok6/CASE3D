//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

int loc;
double** WAVE_IN(int NNODE, double** GCOORD, double* T, int TIME, double** PIN, double DT, double PPEAK, double TAU, double XC, double YC, double ZC, double XO, double YO, double ZO) {
	//In this code, we've only coded the spherical wave case since you don't have plane wave in real application
	int i, j;
	double R = 0.0; 
	double RO = 0.0; 

	if (TIME == 1) {
		RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
		//the distance between reference point and explosion point
		//However, in spherical wave case, SONODE should be the point that spherical wave reach first when contacting with the free surface.  
		for (j = 0; j < NNODE; j++) {
			R = sqrt(pow((GCOORD[j][2] - ZC), 2) + pow((GCOORD[j][0] - XC), 2) + pow((GCOORD[j][1] - YC), 2));
			if ((T[TIME - 1]) - (R - RO) / C >= 0) { //error prone (please refer to the comment in PWA)
			//if (T[TIME - 1] - (R - RO) / C > 1e-6) {
				PIN[j][0] = PPEAK*(RO / R)*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
				//if (R == RO) {
					//count += 1;
				//}
			}
			else {
				PIN[j][0] = 0.0;
			}
		}
		std::cout << " " << std::endl;
	}
	else { //TIME>1
		RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
		for (j = 0; j < NNODE; j++) {
			R = sqrt(pow((GCOORD[j][2] - ZC), 2) + pow((GCOORD[j][0] - XC), 2) + pow((GCOORD[j][1] - YC), 2));
			if ((T[TIME - 1]) - (R - RO) / C >= 0) {
				//if (T[TIME - 1] - (R - RO) / C > 1e-6) {
				PIN[j][1] = PPEAK*(RO / R)*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
			}
			else {
				PIN[j][1] = 0.0;
			}
		}
	}
	return (PIN);
}