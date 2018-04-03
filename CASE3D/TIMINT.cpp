//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
//TIMINT returns time step based on the CFL limit and the number of time steps needed to reach TTERM

//---------------------OUTPUT VARIABLES---------------------//
struct TIMINTstruct TIMINT(double LMAX) {
	//----------------------BEGIN TIMINT COMP.---------------//
	//CALCULATE TIMESTEP	BASED ON MAX. EIGENVAULE AND CFL NUMBER
	TIMINTstruct t;
	t.DT = (CFLFRAC * 2) / (C*sqrt(LMAX*(1 + 2 * BETA)));
	t.DT = t.DT / dtscale;
	/*
	if (FEM == 1) {
		t.DT = XHE / C / sqrt(1 + 2 * BETA); //this would actually make the FEM solution unstable
	}
	*/
	//CALCULATE THE NUMBER OF TIMESTEPS
	//t.NDT = int(TTERM / t.DT) + 1;
	//t.NDT = round(TTERM / t.DT) + 1;
	if (debug == 1) {
		t.DT = 1e-06;
	}
	t.NDT = floor(TTERM / t.DT) + 1;
	std::cout << "DT is: " << t.DT << std::endl;
	std::cout << "number of time step is: " << t.NDT << std::endl;
	//-------------------------END TIMINT COMP.-----------------//
	return t;
}