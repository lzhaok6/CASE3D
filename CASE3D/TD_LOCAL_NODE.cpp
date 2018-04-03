//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

//the definition of 2D LNA arrary
struct TD_LOCAL_NODEstruct TD_LOCAL_NODE(int n) {
	TD_LOCAL_NODEstruct t;
	int i, j, ct;
	t.LNA = new int*[n + 1];
	for (i = 0; i < n + 1; i++) {
		t.LNA[i] = new int[n + 1];
	}

	for (i = 0; i < n + 1; i++) {
		for (j = 0; j < n + 1; j++) {
			t.LNA[i][j] = 0;
		}
	}

	if (n == 1) {
		t.LNA[0][0] = 1;
		t.LNA[1][0] = 2;
		t.LNA[1][1] = 3;
		t.LNA[0][1] = 4;
	}
	if (n == 2) {
		t.LNA[0][0] = 1;
		t.LNA[2][0] = 2;
		t.LNA[2][2] = 3;
		t.LNA[0][2] = 4;
		t.LNA[1][0] = 5;
		t.LNA[2][1] = 6;
		t.LNA[1][2] = 7;
		t.LNA[0][1] = 8;
		t.LNA[1][1] = 9;
	}
	if (n == 3) {
		t.LNA[0][0] = 1;
		t.LNA[3][0] = 2;
		t.LNA[3][3] = 3;
		t.LNA[0][3] = 4;
		t.LNA[1][0] = 5;
		t.LNA[2][0] = 6;
		t.LNA[3][1] = 7;
		t.LNA[3][2] = 8;
		t.LNA[2][3] = 9;
		t.LNA[1][3] = 10;
		t.LNA[0][2] = 11;
		t.LNA[0][1] = 12;
		t.LNA[1][1] = 13;
		t.LNA[2][1] = 14;
		t.LNA[2][2] = 15;
		t.LNA[1][2] = 16;
	}

	if (n > 3) {
		ct = 0;
		for (i = 0; i < n + 1; i++) {
			for (j = 0; j < n + 1; j++) {
				t.LNA[i][j] = ct + 1;
				ct += 1;
			}
		}
	}

	return t;
}