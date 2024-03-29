//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;

//-----------------OUTPUT VARIABLES---------------------------//


struct LOBATTOstruct LOBATTO(int n) {
	LOBATTOstruct t;
	int i;
	int NLOBATTO;
	NLOBATTO = n - 1;  //TOTALLY n-1 NODES 
	t.Z = new double[n - 1];
	t.WL = new double[n - 1];
	for (i = 0; i < n - 1; i++) {
		t.Z[i] = 0.0;
		t.WL[i] = 0.0;
	}

	//-------------WRITE INTEGRATION POINTS AND WEIGHTS-------------//
	if (NLOBATTO == 1) {
		t.Z[0] = 0.0;

		t.WL[0] = 4.0 / 3.0;
	}
	else if (NLOBATTO == 2) {
		t.Z[0] = -1.0 / (sqrt(5.0));
		t.Z[1] = -t.Z[0];

		t.WL[0] = 5.0 / 6.0;
		t.WL[1] = t.WL[0];
	}

	else if (NLOBATTO == 3) {
		t.Z[0] = -sqrt(3.0 / 7.0);
		t.Z[1] = 0.0;
		t.Z[2] = -t.Z[0];

		t.WL[0] = 49.0 / 90.0;
		t.WL[1] = 32.0 / 45.0;
		t.WL[2] = t.WL[0];
	}
	else if (NLOBATTO == 4) {
		t.Z[0] = -0.76505532392946;
		t.Z[1] = -0.28523151648064;
		t.Z[2] = -t.Z[1];
		t.Z[3] = -t.Z[0];

		t.WL[0] = 0.37847495629785;
		t.WL[1] = 0.55485837703549;
		t.WL[2] = t.WL[1];
		t.WL[3] = t.WL[0];
	}
	else if (NLOBATTO == 5) {
		t.Z[0] = -0.83022389627857;
		t.Z[1] = -0.46884879347071;
		t.Z[2] = 0.0;
		t.Z[3] = -t.Z[1];
		t.Z[4] = -t.Z[0];

		t.WL[0] = 0.27682604736157;
		t.WL[1] = 0.43174538120986;
		t.WL[2] = 0.48761904761905;
		t.WL[3] = t.WL[1];
		t.WL[4] = t.WL[0];
	}
	else if (NLOBATTO == 6) {
		t.Z[0] = -0.87174014850961;
		t.Z[1] = -0.59170018143314;
		t.Z[2] = -0.20929921790248;
		t.Z[3] = -t.Z[2];
		t.Z[4] = -t.Z[1];
		t.Z[5] = -t.Z[0];

		t.WL[0] = 0.21070422714350;
		t.WL[1] = 0.34112269248350;
		t.WL[2] = 0.41245879465870;
		t.WL[3] = t.WL[2];
		t.WL[4] = t.WL[1];
		t.WL[5] = t.WL[0];
	}
	else if (NLOBATTO == 7) {
		t.Z[0] = -0.89975799541146;
		t.Z[1] = -0.67718627951074;
		t.Z[2] = -0.36311746382618;
		t.Z[3] = 0.0;
		t.Z[4] = -t.Z[2];
		t.Z[5] = -t.Z[1];
		t.Z[6] = -t.Z[0];

		t.WL[0] = 0.16549536156081;
		t.WL[1] = 0.27453871250016;
		t.WL[2] = 0.34642851097305;
		t.WL[3] = 0.37151927437642;
		t.WL[4] = t.WL[2];
		t.WL[5] = t.WL[1];
		t.WL[6] = t.WL[0];
	}
	else if (NLOBATTO == 8) {
		t.Z[0] = -0.91953390816646;
		t.Z[1] = -0.73877386510551;
		t.Z[2] = -0.47792494981044;
		t.Z[3] = -0.16527895766639;
		t.Z[4] = -t.Z[3];
		t.Z[5] = -t.Z[2];
		t.Z[6] = -t.Z[1];
		t.Z[7] = -t.Z[0];

		t.WL[0] = 0.13330599085107;
		t.WL[1] = 0.22488934206313;
		t.WL[2] = 0.29204268367968;
		t.WL[3] = 0.32753976118390;
		t.WL[4] = t.WL[3];
		t.WL[5] = t.WL[2];
		t.WL[6] = t.WL[1];
		t.WL[7] = t.WL[0];
	}
	else if (NLOBATTO == 9) {
		t.Z[0] = -0.93400143040806;
		t.Z[1] = -0.78448347366314;
		t.Z[2] = -0.56523532699621;
		t.Z[3] = -0.29575813558694;
		t.Z[4] = 0.0;
		t.Z[5] = -t.Z[3];
		t.Z[6] = -t.Z[2];
		t.Z[7] = -t.Z[1];
		t.Z[8] = -t.Z[0];

		t.WL[0] = 0.10961227326699;
		t.WL[1] = 0.18716988178031;
		t.WL[2] = 0.24804810426403;
		t.WL[3] = 0.28687912477901;
		t.WL[4] = 0.30021759545569;
		t.WL[5] = t.WL[3];
		t.WL[6] = t.WL[2];
		t.WL[7] = t.WL[1];
		t.WL[8] = t.WL[0];
	}
	else if (NLOBATTO == 10) {
		t.Z[0] = -0.94489927222288;
		t.Z[1] = -0.81927932164401;
		t.Z[2] = -0.63287615303186;
		t.Z[3] = -0.39953094096535;
		t.Z[4] = -0.13655293285493;
		t.Z[5] = -t.Z[4];
		t.Z[6] = -t.Z[3];
		t.Z[7] = -t.Z[2];
		t.Z[8] = -t.Z[1];
		t.Z[9] = -t.Z[0];

		t.WL[0] = 0.09168451741320;
		t.WL[1] = 0.15797470556437;
		t.WL[2] = 0.21250841776102;
		t.WL[3] = 0.25127560319920;
		t.WL[4] = 0.27140524091070;
		t.WL[5] = t.WL[4];
		t.WL[6] = t.WL[3];
		t.WL[7] = t.WL[2];
		t.WL[8] = t.WL[1];
		t.WL[9] = t.WL[0];
	}
	else if (NLOBATTO == 11) {
		t.Z[0] = -0.95330984664216;
		t.Z[1] = -0.84634756465187;
		t.Z[2] = -0.68618846908176;
		t.Z[3] = -0.48290982109134;
		t.Z[4] = -0.24928693010624;
		t.Z[5] = 0.0;
		t.Z[6] = -t.Z[4];
		t.Z[7] = -t.Z[3];
		t.Z[8] = -t.Z[2];
		t.Z[9] = -t.Z[1];
		t.Z[10] = -t.Z[0];

		t.WL[0] = 0.07780168674682;
		t.WL[1] = 0.13498192668961;
		t.WL[2] = 0.18364686520355;
		t.WL[3] = 0.22076779356611;
		t.WL[4] = 0.24401579030668;
		t.WL[5] = 0.25193084933345;
		t.WL[6] = t.WL[4];
		t.WL[7] = t.WL[3];
		t.WL[8] = t.WL[2];
		t.WL[9] = t.WL[1];
		t.WL[10] = t.WL[0];
	}
	else if (NLOBATTO == 12) {
		t.Z[0] = -0.95993504526726;
		t.Z[1] = -0.86780105383035;
		t.Z[2] = -0.72886859909133;
		t.Z[3] = -0.55063940292865;
		t.Z[4] = -0.34272401334271;
		t.Z[5] = -0.11633186888370;
		t.Z[6] = -t.Z[5];
		t.Z[7] = -t.Z[4];
		t.Z[8] = -t.Z[3];
		t.Z[9] = -t.Z[2];
		t.Z[10] = -t.Z[1];
		t.Z[11] = -t.Z[0];

		t.WL[0] = 0.06683728449768;
		t.WL[1] = 0.11658665589871;
		t.WL[2] = 0.16002185176295;
		t.WL[3] = 0.19482614937342;
		t.WL[4] = 0.21912625300977;
		t.WL[5] = 0.23161279446846;
		t.WL[6] = t.WL[5];
		t.WL[7] = t.WL[4];
		t.WL[8] = t.WL[3];
		t.WL[9] = t.WL[2];
		t.WL[10] = t.WL[1];
		t.WL[11] = t.WL[0];
	}
	else if (NLOBATTO == 13) {
		t.Z[0] = -0.96524592650384;
		t.Z[1] = -0.88508204422298;
		t.Z[2] = -0.76351968995182;
		t.Z[3] = -0.60625320546985;
		t.Z[4] = -0.42063805471367;
		t.Z[5] = -0.21535395536379;
		t.Z[6] = 0.0;
		t.Z[7] = -t.Z[5];
		t.Z[8] = -t.Z[4];
		t.Z[9] = -t.Z[3];
		t.Z[10] = -t.Z[2];
		t.Z[11] = -t.Z[1];
		t.Z[12] = -t.Z[0];

		t.WL[0] = 0.05802989302860;
		t.WL[1] = 0.10166007032572;
		t.WL[2] = 0.14051169980243;
		t.WL[3] = 0.17278964725360;
		t.WL[4] = 0.19698723596461;
		t.WL[5] = 0.21197358592682;
		t.WL[6] = 0.21704811634882;
		t.WL[7] = t.WL[5];
		t.WL[8] = t.WL[4];
		t.WL[9] = t.WL[3];
		t.WL[10] = t.WL[2];
		t.WL[11] = t.WL[1];
		t.WL[12] = t.WL[0];
	}

	else if (NLOBATTO == 14) {
		t.Z[0] = -0.96956804627022;
		t.Z[1] = -0.89920053309347;
		t.Z[2] = -0.79200829186182;
		t.Z[3] = -0.65238870288249;
		t.Z[4] = -0.48605942188714;
		t.Z[5] = -0.29983046890076;
		t.Z[6] = -0.10132627352195;
		t.Z[7] = -t.Z[6];
		t.Z[8] = -t.Z[5];
		t.Z[9] = -t.Z[4];
		t.Z[10] = -t.Z[3];
		t.Z[11] = -t.Z[2];
		t.Z[12] = -t.Z[1];
		t.Z[13] = -t.Z[0];

		t.WL[0] = 0.05085036100592;
		t.WL[1] = 0.08939369732593;
		t.WL[2] = 0.12425538213251;
		t.WL[3] = 0.15402698080716;
		t.WL[4] = 0.17749191339170;
		t.WL[5] = 0.19369002382520;
		t.WL[6] = 0.20195830817823;
		t.WL[7] = t.WL[6];
		t.WL[8] = t.WL[5];
		t.WL[9] = t.WL[4];
		t.WL[10] = t.WL[3];
		t.WL[11] = t.WL[2];
		t.WL[12] = t.WL[1];
		t.WL[13] = t.WL[0];
	}
	else if (NLOBATTO == 15) {
		t.Z[0] = -0.97313217663142;
		t.Z[1] = -0.91087999591557;
		t.Z[2] = -0.81569625122177;
		t.Z[3] = -0.69102898062768;
		t.Z[4] = -0.54138539933010;
		t.Z[5] = -0.37217443356548;
		t.Z[6] = -0.18951197351832;
		t.Z[7] = 0.0;
		t.Z[8] = -t.Z[6];
		t.Z[9] = -t.Z[5];
		t.Z[10] = -t.Z[4];
		t.Z[11] = -t.Z[3];
		t.Z[12] = -t.Z[2];
		t.Z[13] = -t.Z[1];
		t.Z[14] = -t.Z[0];

		t.WL[0] = 0.04492194054325;
		t.WL[1] = 0.07919827050369;
		t.WL[2] = 0.11059290900703;
		t.WL[3] = 0.13798774620193;
		t.WL[4] = 0.16039466199762;
		t.WL[5] = 0.17700425351566;
		t.WL[6] = 0.18721633967762;
		t.WL[7] = 0.19066187475347;
		t.WL[8] = t.WL[6];
		t.WL[9] = t.WL[5];
		t.WL[10] = t.WL[4];
		t.WL[11] = t.WL[3];
		t.WL[12] = t.WL[2];
		t.WL[13] = t.WL[1];
		t.WL[14] = t.WL[0];
	}
	return t;
}