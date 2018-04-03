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
double** WAVE_IN(int NNODE, double** GCOORD, double* T, int TIME, double** PIN, int *NRBA, int NRBNODE, double*timer, double*ampt, double DT, double PPEAK, double TAU, double XC, double YC, double ZC, double XO, double YO, double ZO) {
	int i, j;
	double R = 0.0; 
	double RO = 0.0; 
	if (WAVE == 1) {
		if (Abaquswaveform == 1) {
			//currently, the abaqus wave form can only be added on the NRB
			if (TIME == 1) {
				for (i = 0; i < NRBNODE; i++) {
					PIN[NRBA[i] - 1][0] = PPEAK*ampt[0];
				}
				loc = 1;
			}
			else {  //TIME>1
				if (T[TIME - 1] < timer[loc]) {
					for (i = 0; i < NRBNODE; i++) {
						PIN[NRBA[i] - 1][1] = PPEAK*ampt[loc - 1] + PPEAK*(ampt[loc] - ampt[loc - 1])*((T[TIME - 1] - timer[loc - 1]) / (timer[loc] - timer[loc - 1]));
					}
				}
				else {
					if (loc < 115) {
						loc += 1;
						for (i = 0; i < NRBNODE; i++) {
							PIN[NRBA[i] - 1][1] = PPEAK*ampt[loc - 1] + PPEAK*(ampt[loc] - ampt[loc - 1])*((T[TIME - 1] - timer[loc - 1]) / (timer[loc] - timer[loc - 1]));
						}
						//the last pair is timer[114] and timer[115]
					}
					else { //if loc==115
						for (i = 0; i < NRBNODE; i++) {
							PIN[NRBA[i] - 1][1] = PPEAK*ampt[loc] + PPEAK*(0 - ampt[loc])*((T[TIME - 1] - timer[loc]) / (TTERM - timer[loc]));
						}
					} //after the last ampt definition, linear
				}
			}
		}
		else { //WAVE==1 and if the exponential decay sharp wave form is used
			if (TIME == 1) {
				for (j = 0; j < NNODE; j++) {
					//if ((T[TIME - 1] + (abs(GCOORD[j][1]) - xo) / C) >= 0) {
					if ((T[TIME - 1] + ((-xo) - GCOORD[j][1]) / C) >= 0) {
						//PIN[j][0] = PPEAK*exp(-(T[TIME - 1] + ((abs(GCOORD[j][1]) - xo) / C)) / TAU);
						PIN[j][0] = PPEAK*exp(-(T[TIME - 1] + (((-xo) - GCOORD[j][1]) / C)) / TAU);
					}
					else {
						PIN[j][0] = 0.0;
					}
				}
			}
			else { //TIME>0
				for (j = 0; j < NNODE; j++) {
					//if ((T[TIME - 1] + (abs(GCOORD[j][1]) - xo) / C) >= 0) {
					if ((T[TIME - 1] + ((-xo) - GCOORD[j][1]) / C) >= 0) {
						//PIN[j][1] = PPEAK*exp(-(T[TIME - 1] + ((abs(GCOORD[j][1]) - xo) / C)) / TAU);
						PIN[j][1] = PPEAK*exp(-(T[TIME - 1] + (((-xo) - GCOORD[j][1]) / C)) / TAU);
					}
					else {
						PIN[j][1] = 0.0;
					}
				}
			}
		}
	}

	//the WAVE==2 case (spherical wave)
	else {
		//int count = 0; 
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
	}
	return (PIN);
}