//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;

//LN[][0] is the value of Legendre polynomial and LN[][1] is the value of the derivative of Legendre polynomial

struct LEGENDREstruct LEGENDRE(int Nq, int n) {
	LEGENDREstruct t;
	int i;
	t.LN = new double*[2];
	for (i = 0; i < 2; i++) {
		t.LN[i] = new double[Nq + 1];
	}

	if (n == 1) {
		if (Nq == n) {
			t.LN[0][0] = -1.0;
			t.LN[0][1] = 1.0;

			t.LN[1][0] = 1.0;
			t.LN[1][1] = t.LN[1][0];
		}
		else if (Nq > n) {
			t.LN[0][0] = -1.0;
			t.LN[0][Nq] = 1.0;
			for (i = 1; i < Nq; i++) {
				t.LN[0][i] = 0.0;
			}
			for (i = 0; i < Nq + 1; i++) {
				t.LN[1][i] = 1.0;
			}
		}
		/*
		else if (Nq == n + 1) { //If one more Lobatto point is inserted
			t.LN[0][0] = -1.0;
			t.LN[0][1] = 0.0;
			t.LN[0][2] = 1.0;

			t.LN[1][0] = 1.0;
			t.LN[1][1] = 1.0; 
			t.LN[1][2] = 1.0;
		}
		else if (Nq == n + 2) {
			t.LN[0][0] = -1.0;
			t.LN[0][1] = 0.0;
			t.LN[0][2] = 0.0;
			t.LN[0][3] = 1.0;

			t.LN[1][0] = 1.0;
			t.LN[1][1] = 1.0;
			t.LN[1][2] = 1.0;
			t.LN[1][3] = 1.0;
		}
		*/
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 2) {
		if (Nq == n) {
			t.LN[0][0] = 1.0;
			t.LN[0][1] = -0.5;
			t.LN[0][2] = t.LN[0][0];

			t.LN[1][0] = -3.0;
			t.LN[1][1] = 0.0;
			t.LN[1][2] = -t.LN[1][0];
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = 1.000000000000000;
			t.LN[0][1] = -0.200000000000000;
			t.LN[0][2] = -0.200000000000000;
			t.LN[0][3] = 1.000000000000000;

			t.LN[1][0] = -3.000000000000000;
			t.LN[1][1] = -1.341640786499874;
			t.LN[1][2] = 1.341640786499874;
			t.LN[1][3] = 3.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 3) {
		if (Nq == n) {
			t.LN[0][0] = -1.0;
			t.LN[0][1] = 0.44721359549996;
			t.LN[0][2] = -t.LN[0][1];
			t.LN[0][3] = -t.LN[0][0];

			t.LN[1][0] = 6.0;
			t.LN[1][1] = 0.0;
			t.LN[1][2] = 0.0;
			t.LN[1][3] = 6.0;
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = -1.000000000000000;
			t.LN[0][1] = 0.280565858874847;
			t.LN[0][2] = 0.0;
			t.LN[0][3] = -0.280565858874847;
			t.LN[0][4] = 1.000000000000000;

			t.LN[1][0] = 6.000000000000000;
			t.LN[1][1] = 1.714285714285714;
			t.LN[1][2] = -1.500000000000000;
			t.LN[1][3] = 1.714285714285714;
			t.LN[1][4] = 6.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 4) {
		if (Nq == n) {
			t.LN[0][0] = 1.0;
			t.LN[0][1] = -0.42857142857143;
			t.LN[0][2] = 0.375;
			t.LN[0][3] = t.LN[0][1];
			t.LN[0][4] = t.LN[0][0];

			t.LN[1][0] = -10.0;
			t.LN[1][1] = 0.0;
			t.LN[1][2] = 0.0;
			t.LN[1][3] = 0.0;
			t.LN[1][4] = -t.LN[1][0];
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = 1.000000000000000;
			t.LN[0][1] = -0.321091373894025;
			t.LN[0][2] = 0.098869151671802;
			t.LN[0][3] = 0.098869151671802;
			t.LN[0][4] = -0.321091373894025;
			t.LN[0][5] = 1.000000000000000;

			t.LN[1][0] = -10.000000000000000;
			t.LN[1][1] = -2.098484670656326;
			t.LN[1][2] = 1.733138625277072;
			t.LN[1][3] = -1.733138625277072;
			t.LN[1][4] = 2.098484670656326;
			t.LN[1][5] = 10.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 5) {
		if (Nq == n) {
			t.LN[0][0] = -1.0;
			t.LN[0][1] = 0.41969693413129;
			t.LN[0][2] = -0.34662772505542;
			t.LN[0][3] = -t.LN[0][2];
			t.LN[0][4] = -t.LN[0][1];
			t.LN[0][5] = -t.LN[0][0];

			t.LN[1][0] = 15.0;
			t.LN[1][1] = -1.421085471520200E-13;
			t.LN[1][2] = 5.773159728050814E-14;
			t.LN[1][3] = t.LN[1][2];
			t.LN[1][4] = t.LN[1][1];
			t.LN[1][5] = t.LN[1][0];
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = -1.000000000000000;
			t.LN[0][1] = 0.344335743198449;
			t.LN[0][2] = -0.155707418768582;
			t.LN[0][3] = 0.0;
			t.LN[0][4] = 0.155707418768582;
			t.LN[0][5] = -0.344335743198449;
			t.LN[0][6] = 1.000000000000000;

			t.LN[1][0] = 15.000000000000000;
			t.LN[1][1] = 2.488502762268895;
			t.LN[1][2] = -1.992634993673677;
			t.LN[1][3] = 1.875000000000000;
			t.LN[1][4] = -1.992634993673677;
			t.LN[1][5] = 2.488502762268895;
			t.LN[1][6] = 15.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 6) {
		if (Nq == n) {
			t.LN[0][0] = 1.0;
			t.LN[0][1] = -0.41475046037813;
			t.LN[0][2] = 0.33210583227895;
			t.LN[0][3] = -0.3125;
			t.LN[0][4] = t.LN[0][2];
			t.LN[0][5] = t.LN[0][1];
			t.LN[0][6] = t.LN[0][0];

			t.LN[1][0] = -21.0;
			t.LN[1][1] = -1.652011860642233E-13;
			t.LN[1][2] = -7.549516567451065E-14;
			t.LN[1][3] = 0.0;
			t.LN[1][4] = -t.LN[1][2];
			t.LN[1][5] = -t.LN[1][1];
			t.LN[1][6] = -t.LN[1][0];
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = 1.000000000000000;
			t.LN[0][1] = -0.358898305790205;
			t.LN[0][2] = 0.191455294719785;
			t.LN[0][3] = -0.061588312635442;
			t.LN[0][4] = -0.061588312635442;
			t.LN[0][5] = 0.191455294719785;
			t.LN[0][6] = -0.358898305790205;
			t.LN[0][7] = 1.000000000000000;

			t.LN[1][0] = -21.000000000000000;
			t.LN[1][1] = -2.881923179547201;
			t.LN[1][2] = 2.264976596411382;
			t.LN[1][3] = -2.059817484119737;
			t.LN[1][4] = 2.059817484119737;
			t.LN[1][5] = -2.264976596411382;
			t.LN[1][6] = 2.881923179547201;
			t.LN[1][7] = 21.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 7) {
		if (Nq == n) {
			t.LN[0][0] = -1.0;
			t.LN[0][1] = 0.41170331136384;
			t.LN[0][2] = -0.32356808520163;
			t.LN[0][3] = 0.29425964058853;
			t.LN[0][4] = -t.LN[0][3];
			t.LN[0][5] = -t.LN[0][2];
			t.LN[0][6] = -t.LN[0][1];
			t.LN[0][7] = -t.LN[0][0];

			t.LN[1][0] = 28.0;
			t.LN[1][1] = 3.197442310920451E-13;
			t.LN[1][2] = 6.750155989720952E-14;
			t.LN[1][3] = 1.953992523340276E-14;
			t.LN[1][4] = t.LN[1][3];
			t.LN[1][5] = t.LN[1][2];
			t.LN[1][6] = t.LN[1][1];
			t.LN[1][7] = t.LN[1][0];
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = -1.000000000000000;
			t.LN[0][1] = 0.368622254677538;
			t.LN[0][2] = -0.215404664516627;
			t.LN[0][3] = 0.102822716279203;
			t.LN[0][4] = 0;
			t.LN[0][5] = -0.102822716279203;
			t.LN[0][6] = 0.215404664516627;
			t.LN[0][7] = -0.368622254677538;
			t.LN[0][8] = 1.000000000000000;

			t.LN[1][0] = 28.000000000000000;
			t.LN[1][1] = 3.277523570181452;
			t.LN[1][2] = -2.544702053588752;
			t.LN[1][3] = 2.265332329561071;
			t.LN[1][4] = -2.187500000000000;
			t.LN[1][5] = 2.265332329561071;
			t.LN[1][6] = -2.544702053588752;
			t.LN[1][7] = 3.277523570181452;
			t.LN[1][8] = 28.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 8) {
		if (Nq == n) {
			t.LN[0][0] = 1.0;
			t.LN[0][1] = -0.40969044627268;
			t.LN[0][2] = 0.31808775669859;
			t.LN[0][3] = -0.28316654119513;
			t.LN[0][4] = 0.27343750000000;
			t.LN[0][5] = t.LN[0][3];
			t.LN[0][6] = t.LN[0][2];
			t.LN[0][7] = t.LN[0][1];
			t.LN[0][8] = t.LN[0][0];

			t.LN[1][0] = -36.0;
			t.LN[1][1] = 2.131628207280301E-14;
			t.LN[1][2] = 8.526512829121202E-14;
			t.LN[1][3] = -4.263256414560601E-14;
			t.LN[1][4] = 0.0;
			t.LN[1][5] = -t.LN[1][3];
			t.LN[1][6] = -t.LN[1][2];
			t.LN[1][7] = -t.LN[1][1];
			t.LN[1][8] = -t.LN[1][0];
		}
		else if (Nq == n + 1) {
			t.LN[0][0] = 1.000000000000000;
			t.LN[0][1] = -0.375436643142831;
			t.LN[0][2] = 0.232231340794229;
			t.LN[0][3] = -0.131834865664879;
			t.LN[0][4] = 0.043050608501474;
			t.LN[0][5] = 0.043050608501474;
			t.LN[0][6] = -0.131834865664879;
			t.LN[0][7] = 0.232231340794229;
			t.LN[0][8] = -0.375436643142831;
			t.LN[0][9] = 1.000000000000000;

			t.LN[1][0] = -36.000000000000000;
			t.LN[1][1] = -3.674611407232717;
			t.LN[1][2] = 2.829122910093557;
			t.LN[1][3] = -2.482636220298543;
			t.LN[1][4] = 2.344251694129179;
			t.LN[1][5] = -2.344251694129179;
			t.LN[1][6] = 2.482636220298543;
			t.LN[1][7] = -2.829122910093557;
			t.LN[1][8] = 3.674611407232717;
			t.LN[1][9] = 36.000000000000000;
		}
		else {
			std::cout << "LN is not defined" << std::endl;
		}
	}
	else if (n == 16) {
		t.LN[0][0] = 1.0;
		t.LN[0][1] = -0.00404577094794E2;
		t.LN[0][2] = 0.00304700170288E2;
		t.LN[0][3] = -0.00257849862958E2;
		t.LN[0][4] = 0.00230839585175E2;
		t.LN[0][5] = -0.00214109329017E2;
		t.LN[0][6] = 0.00203816197260E2;
		t.LN[0][7] = -0.00198179473745E2;
		t.LN[0][8] = 0.00196380615234E2;
		t.LN[0][9] = t.LN[0][7];
		t.LN[0][10] = t.LN[0][6];
		t.LN[0][11] = t.LN[0][5];
		t.LN[0][12] = t.LN[0][4];
		t.LN[0][13] = t.LN[0][3];
		t.LN[0][14] = t.LN[0][2];
		t.LN[0][15] = t.LN[0][1];
		t.LN[0][16] = t.LN[0][0];

		t.LN[1][0] = -1.36E2;
		t.LN[1][1] = -1.592326270838385E-11;
		t.LN[1][2] = -1.634248292248230E-13;
		t.LN[1][3] = -3.346656285430072E-12;
		t.LN[1][4] = 1.136868377216160E-13;
		t.LN[1][5] = -4.334310688136611E-13;
		t.LN[1][6] = 1.847411112976261E-13;
		t.LN[1][7] = -1.438849039914203E-13;
		t.LN[1][8] = 0.0;
		t.LN[1][9] = -t.LN[1][7];
		t.LN[1][10] = -t.LN[1][6];
		t.LN[1][11] = -t.LN[1][5];
		t.LN[1][12] = -t.LN[1][4];
		t.LN[1][13] = -t.LN[1][3];
		t.LN[1][14] = -t.LN[1][2];
		t.LN[1][15] = -t.LN[1][1];
		t.LN[1][16] = -t.LN[1][0];
	}
	return t;
}