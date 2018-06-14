#include "stdafx.h"
#include "Header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

struct LOCAL_NODEstruct LOCAL_NODE(int n) {
	int CN;
	int IC;  //loop counter
	int i, j, k, m;
	LOCAL_NODEstruct t;
	t.LNA = new int**[n + 1];
	for (i = 0; i < n + 1; i++) {
		t.LNA[i] = new int*[n + 1];
		for (j = 0; j < n + 1; j++) {
			t.LNA[i][j] = new int[n + 1];
		}
	}

	for (i = 0; i < n + 1; i++) {
		for (j = 0; j < n + 1; j++) {
			for (k = 0; k < n + 1; k++) {
				t.LNA[i][j][k] = 0;
			}
		}
	}

	//if (n == 1) {
	t.LNA[0][0][0] = 1;
	t.LNA[n][0][0] = 2;
	t.LNA[n][n][0] = 3;
	t.LNA[0][n][0] = 4;
	t.LNA[0][0][n] = 5;
	t.LNA[n][0][n] = 6;
	t.LNA[n][n][n] = 7;
	t.LNA[0][n][n] = 8;
	//}
	//import mesh from outside
	if (n > 1) {
		int u, v, w;
		int NNODE;
		int *ELEMENT_POINT;
		double **FGCOORD;
		std::string filename = "LNA_info_N=" + std::to_string(N) + ".msh";
		std::ifstream infile(filename);   //LNA_info_N=x.txt is for xth order LNA
		if (!infile) {
			std::cerr << "can not open the LNA information file" << std::endl;
			exit(1);
		}
		//read the file 
		int nodestart = 0;
		int nodeend = 0;
		int ele_line = 0;
		std::vector<std::vector<std::string>> output;
		int ct = -1;
		std::string csvLine;
		while (getline(infile, csvLine))
		{
			ct = ct + 1; //the current line number (starting from 0)
			std::istringstream csvStream(csvLine);
			std::vector<std::string> csvColumn;
			std::string csvElement;
			while (getline(csvStream, csvElement, ' '))
			{
				csvColumn.push_back(csvElement);
			}
			output.push_back(csvColumn);

			if (csvColumn[0] == "$Nodes") {
				nodestart = ct + 2; //the node starts from the next line
			}
			if (csvColumn[0] == "$EndNodes") {
				nodeend = ct - 1; //the node starts from the next line
			}
			if (csvColumn[0] == "$Elements") {
				ele_line = ct + 2; //the line corresponding to element connectivity
			}
		}

		NNODE = stoi(output[nodestart - 1][0]);
		ELEMENT_POINT = new int[NINT*NINT*NINT]; //Connectivity matrix
		for (i = 0; i < NINT*NINT*NINT; i++) {
			ELEMENT_POINT[i] = 0;
		}
		for (i = 0; i < NINT*NINT*NINT; i++) {
			ELEMENT_POINT[i] = stoi(output[ele_line][6 + i]);
		}

		FGCOORD = new double*[NNODE];
		for (i = 0; i < NNODE; i++) {
			FGCOORD[i] = new double[3];
		}

		//In the current LNA template file, the origin is (-1,-1,0). It is supposed to start at (0,0,0). We shift the coordinate here but need to 
		//be alerted that we might need to change it in the future.
		for (i = 0; i < NNODE; i++) {
			for (j = 0; j < 3; j++) {
				if (j == 0) {
					FGCOORD[i][j] = stod(output[nodestart + i][1 + j]);
				}
				if (j == 1) {
					FGCOORD[i][j] = stod(output[nodestart + i][1 + j]);
				}
				if (j == 2) {
					FGCOORD[i][j] = stod(output[nodestart + i][1 + j]);
				}
			}
		}

		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < 3; j++) {
				t.L[i][j] = 0;
			}
		}

		double elementsize = 1.0;
		for (i = 0; i < NINT*NINT*NINT; i++) {
			u = round(FGCOORD[ELEMENT_POINT[i] - 1][0] / (elementsize / N));
			v = round(FGCOORD[ELEMENT_POINT[i] - 1][1] / (elementsize / N));
			w = round(FGCOORD[ELEMENT_POINT[i] - 1][2] / (elementsize / N));
			t.L[i][0] = u; t.L[i][1] = v; t.L[i][2] = w;
			t.LNA[u][v][w] = i + 1;
		}
	}
	std::cout << " " << std::endl;
	return t;
}