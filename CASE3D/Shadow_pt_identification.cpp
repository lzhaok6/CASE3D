#include "header.h"
#include "data.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>
#include <sstream>
#include <cstring>
#include <algorithm> 


//This code identifies shadow points in the fluid domain. Only applicable to SFM. 
//This function loop through the fluid domain points. If the line segment connecting the point and charge center intersects. 
//one of the fluid wetted surface elements, then we can say that the fluid node is in the shadow.
//Check out the post-processing folder shadow1.jpg and shadow2.jpg on how to determine if line and element intersect. 
std::vector<int> shadow_pt_identification(double XC, double YC, double ZC, int NNODE, double** GCOORD) {
	extern STRU_WET_SURF ss[ssnumber];
	extern OWETSURF ol[owsfnumber];
	//First determine the potential points in the shadow (the assumption is that the center line is at z=0)
	//determine the Xmax and Xmin
	std::vector<double> xcoord; 
	std::vector<double> ycoord;
	std::vector<double> zcoord;
	for (int z = 0; z < ssnumber; z++) {
		for (int i = 0; i < ss[z].Node_stru; i++) {
			xcoord.push_back(ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][0]);
			ycoord.push_back(ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][1]);
			zcoord.push_back(ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][2]);
		}
	}
	std::sort(xcoord.begin(), xcoord.end());
	std::sort(ycoord.begin(), ycoord.end());
	std::sort(zcoord.begin(), zcoord.end());
	double Xmin = xcoord.front();
	double Xmax = xcoord.back();
	double Ymin = ycoord.front();
	double Ymax = ycoord.back();
	double Zmin = zcoord.front();
	double Zmax = zcoord.back();
	//The shadow points should be confined within the region of [Xmin~Xmax, Ymin~Ymax, 0~Zmax/Zmin~0 depending on which side is the charge]
	if (ZC > 0) { // on the other side of the charge
		Zmax = 0.0; 
	}
	else {
		Zmin = 0.0; 
	}
	std::vector<int> shadow_pt_raw; 
	//loop through the fluid nodes to find the tentative shadow region
	for (int i = 0; i < NNODE; i++) {
		if (ZC > 0) {
			if (GCOORD[i][1] >= Ymin && GCOORD[i][1] <= Ymax && GCOORD[i][2] <= 0.0) {
				shadow_pt_raw.push_back(i + 1);
			}
		}
		else {
			if (GCOORD[i][1] >= Ymin && GCOORD[i][1] <= Ymax && GCOORD[i][2] >= 0.0) {
				shadow_pt_raw.push_back(i + 1);
			}
		}
	}
	
	//Begin to search for the real shadow points by determining if the fluid and charge center point intersects
	std::vector<int> shadow_pts;
	//Now we search for the points within the box
	double V0[3];
	double P0[3];
	double P1[3];
	P0[0] = XC; P0[1] = YC; P0[2] = ZC;
	double pt_inter[3]; //point of intersection
	int inelement;
	double xi, eta;
	int intersection; 
	double elexmin, elexmax, eleymin, eleymax, elezmin, elezmax; 
	double elex[4], eley[4], elez[4];
	int elenode2D; 
	if (element_type == 0) {
		elenode2D = 4; 
	}
	else if (element_type == 1) {
		elenode2D = 3; 
	}
	//loop through the fluid points
	for (int i = 0; i < shadow_pt_raw.size(); i++) {
		inelement = 0; //initialize assuming no intersection
		P1[0] = GCOORD[shadow_pt_raw[i] - 1][0]; P1[1] = GCOORD[shadow_pt_raw[i] - 1][1]; P1[2] = GCOORD[shadow_pt_raw[i] - 1][2];
		for (int z = 0; z < ssnumber; z++) {
			for (int j = 0; j < ol[z].FSNEL; j++) {
				//loop through structural elements and see if the line intersects with any of the element. 
				//If they intersect, then the point is in shadow, jump out of the loop of structural element. 
				V0[0] = GCOORD[ol[z].IEN_algo2[0][j] - 1][0];
				V0[1] = GCOORD[ol[z].IEN_algo2[0][j] - 1][1];
				V0[2] = GCOORD[ol[z].IEN_algo2[0][j] - 1][2];
				double s = (ol[z].norm[j][0] * (V0[0] - P0[0]) + ol[z].norm[j][1] * (V0[1] - P0[1]) + ol[z].norm[j][2] * (V0[2] - P0[2])) /
					(ol[z].norm[j][0] * (P1[0] - P0[0]) + ol[z].norm[j][1] * (P1[1] - P0[1]) + ol[z].norm[j][2] * (P1[2] - P0[2]));
				if (s >=0 && s < 0.99) { 
					//the line has intersection with the element plane. We can move on to test if the point is within the element.
					//the reason to use s<0.99 is to make sure that the point inside its corresponding node would not be considered. (s==1 means the intersection is at the end of line segment)
					intersection = 1;
				}
				else { //If there is even no intersection, we can move on to find another element. 
					intersection = 0;
				}
				if (intersection == 1) {
					//Get the point of intersection (coordinate)
					pt_inter[0] = P0[0] + s*(P1[0] - P0[0]); pt_inter[1] = P0[1] + s*(P1[1] - P0[1]); pt_inter[2] = P0[2] + s*(P1[2] - P0[2]);
					//next determine if the point is inside the element
					//since we don't have to determine the local coordinate of the projected point, we can use a simpler method to determine if the point is in the element
					//First, we find the box bounds the 4-node element [Xmin~Xmax, Ymin~Ymax, Zmin~Zmax]. If the point is in the bounding box, then it means the projected node is within the element. 

					for (int k = 0; k < elenode2D; k++) {
						elex[k] = GCOORD[ol[z].IEN_algo2[k][j] - 1][0];
						eley[k] = GCOORD[ol[z].IEN_algo2[k][j] - 1][1];
						elez[k] = GCOORD[ol[z].IEN_algo2[k][j] - 1][2];
					}
					if (elex[0] > elex[1]) {
						elexmax = elex[0]; elexmin = elex[1];
					}
					else {
						elexmin = elex[0]; elexmax = elex[1];
					}
					if (eley[0] > eley[1]) {
						eleymax = eley[0]; eleymin = eley[1];
					}
					else {
						eleymin = eley[0]; eleymax = eley[1];
					}
					if (elez[0] > elez[1]) {
						elezmax = elez[0]; elezmin = elez[1];
					}
					else {
						elezmin = elez[0]; elezmax = elez[1];
					}
					for (int k = 2; k < elenode2D; k++) {
						if (elex[k] > elexmax) {
							elexmax = elex[k];
						}
						else if (elex[k] < elexmin) {
							elexmin = elex[k];
						}
						if (eley[k] > eleymax) {
							eleymax = eley[k];
						}
						else if (eley[k] < eleymin) {
							eleymin = eley[k];
						}
						if (elez[k] > elezmax) {
							elezmax = elez[k];
						}
						else if (elez[k] < elezmin) {
							elezmin = elez[k];
						}
					}
					if (pt_inter[0] >= elexmin && pt_inter[0] <= elexmax
						&& pt_inter[1] >= eleymin && pt_inter[1] <= eleymax
						&& pt_inter[2] >= elezmin && pt_inter[2] <= elezmax) {
						inelement = 1; //the intersection point is indeed in the element
						goto endsearch;
					}
				}
			}
		}
	endsearch:;
		if (inelement == 1) {
			//The point is not in the shadow region, record the point to be used in the future. 
			shadow_pts.push_back(shadow_pt_raw[i]);
		}
		
	}
	endloop:;

	std::cout << " " << std::endl; 

	return shadow_pts;
}