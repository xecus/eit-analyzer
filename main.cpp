//std
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <ncurses.h>
#include <cctype>
#include <algorithm>

//eigen
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include "largemat.hpp"
#include "gauss.hpp"
#include "domain.hpp"
#include "electrode.hpp"
#include "element.hpp"
#include "fem.hpp"
#include "inverse.hpp"
#include "util.hpp"

/* for ncurses */
void wCenterTitle(WINDOW *pwin, const char * title){
	int x, maxy, maxx, stringsize;
	getmaxyx(pwin, maxy, maxx);
	stringsize = 4 + strlen(title);
	x = (maxx - stringsize)/2;
	mvwaddch(pwin, 0, x, ACS_RTEE);
	waddch(pwin, ' ');
	waddstr(pwin, title);
	waddch(pwin, ' ');
	waddch(pwin, ACS_LTEE);
}
void wclrscr(WINDOW * pwin){
	int y, x, maxy, maxx;
	getmaxyx(pwin, maxy, maxx);
	for(y=0; y < maxy; y++)
		for(x=0; x < maxx; x++)
			mvwaddch(pwin, y, x, ' ');
}


int main1(int argc,char *argv[]){

	//Welcome Message
	std::cout << "[EIT Reconstruction Program]" << std::endl;

	/*
	ESet *eset = new ESet(
		//(new AnalysisDomain())->Set(2,(int[]){8,8,0},(double[]){0.15,0.15,0.08},(int[]){0,0,1})
		(new AnalysisDomain())->Test2D_Exp1()
	);
	eset->Print();
	*/

	InverseProblem *p = new InverseProblem(
		new Fem(
			new ESet(
				(new AnalysisDomain())->Test2D_Exp1(),
				"uho.csv",594
			)
		)
	);

	p->Init();

	{
		Eigen::VectorXd SIGMA = Eigen::VectorXd::Zero(p->femp->SumOfElement);
		for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA(i) = 1.0;
		Eigen::VectorXd tmp = p->femp->CalcZ( SIGMA );

		std::ofstream ofs("Z.csv");
		for(int i=0;i<(p->femp->esp->NumOfES);i++){
			ofs << tmp(i) << std::endl;
		}
		ofs.close();

	}

	return 0;
}


int main3(int argc,char *argv[]){

	//Welcome Message
	std::cout << "[EIT Reconstruction Program]" << std::endl;

	InverseProblem *p = new InverseProblem(new Fem(new ESet(
		//(new AnalysisDomain())->Set(2,(int[]){8,8,0},(double[]){0.15,0.15,0.08},(int[]){0,0,1})
		(new AnalysisDomain())->Test2D_Exp1()
	)));

	std::cout << "[Init()]" << std::endl;
	p->Init();


	int SampleNum = 64*2 + 1;
	int SampleCnt = 0;

	Eigen::MatrixXd SIGMA = Eigen::MatrixXd::Zero( SampleNum , p->femp->SumOfElement);
	Eigen::MatrixXd ZZZ = Eigen::MatrixXd::Zero(  SampleNum , p->femp->esp->NumOfES );

	
	for(int y=0;y<(p->femp->esp->adp->NumOfEleY);y++){
	for(int x=0;x<(p->femp->esp->adp->NumOfEleX);x++){
		std::cout << "A[" << x << "," << y << "]" << std::endl;
		for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA( SampleCnt , i ) = 1.0;
		SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(x,y) ) = 2.0;
		Eigen::VectorXd tmp = SIGMA.row(SampleCnt);
		ZZZ.row(SampleCnt) = p->femp->CalcZ( tmp );
		SampleCnt++;
	}
	}
	
	/*
	for(int y=0;y<(p->femp->esp->adp->NumOfEleY);y++){
	for(int x=0;x<(p->femp->esp->adp->NumOfEleX);x++){
		std::cout << "B[" << x << "," << y << "]" << std::endl;
		for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA( SampleCnt , i ) = 1.0;
		SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(x,y) ) = 5.0;
		Eigen::VectorXd tmp = SIGMA.row(SampleCnt);
		ZZZ.row(SampleCnt) = p->femp->CalcZ( tmp );
		SampleCnt++;
	}
	}
	*/

	//Flat
	{
	for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA( SampleCnt , i ) = 1.0;
	Eigen::VectorXd tmp = SIGMA.row(SampleCnt);
	ZZZ.row(SampleCnt) = p->femp->CalcZ( tmp );
	SampleCnt++;
	}

	//1
	{
	for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA( SampleCnt , i ) = 1.0;
	SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(1,1) ) = 2.0;
	SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(6,6) ) = 5.0;
	Eigen::VectorXd tmp = SIGMA.row(SampleCnt);
	ZZZ.row(SampleCnt) = p->femp->CalcZ( tmp );
	SampleCnt++;
	}

	//1
	{
	for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA( SampleCnt , i ) = 1.0;
	SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(1,1) ) = 2.0;
	SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(6,1) ) = 5.0;
	Eigen::VectorXd tmp = SIGMA.row(SampleCnt);
	ZZZ.row(SampleCnt) = p->femp->CalcZ( tmp );
	SampleCnt++;
	}

	//1
	{
	for(int i=0;i<(p->femp->SumOfElement);i++) SIGMA( SampleCnt , i ) = 1.0;
	SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(1,6) ) = 2.0;
	SIGMA( SampleCnt , p->femp->esp->adp->ConvertElement(6,6) ) = 5.0;
	Eigen::VectorXd tmp = SIGMA.row(SampleCnt);
	ZZZ.row(SampleCnt) = p->femp->CalcZ( tmp );
	SampleCnt++;
	}


	std::ofstream ofs1("./tmp/S.csv");
	std::ofstream ofs2("./tmp/Z.csv");
	for(int i=0;i<SampleCnt;i++){
		for(int j=0;j<(p->femp->SumOfElement);j++){
			ofs1 << SIGMA( i , j ) ;
			if( j !=(p->femp->SumOfElement-1) ) ofs1 << ","; else ofs1 << std::endl;
		}
		for(int j=0;j<(p->femp->esp->NumOfES);j++){
			ofs2 << ZZZ( i , j ) ;
			if( j !=(p->femp->esp->NumOfES-1) ) ofs2 << ","; else ofs2 << std::endl;
		}
	}
	ofs1.close();
	ofs2.close();

	return 0;
}


int main2(int argc,char *argv[]){

	//Welcome Message
	std::cout << "[EIT Reconstruction Program]" << std::endl;

	InverseProblem *p = new InverseProblem(new Fem(new ESet(
		//(new AnalysisDomain())->Set(2,(int[]){8,8,0},(double[]){0.15,0.15,0.08},(int[]){0,0,1})
		(new AnalysisDomain())->Test2D_Exp1()
	)));

	std::cout << "[Init()]" << std::endl;
	p->Init();

	std::cout << "[LoadPhantomZ()]" << std::endl;
	p->LoadPhantomZ();

	//std::cout << p->MakeKanzanFilter() << std::endl;

	for(int cnt=0;;cnt++){
		//Show
		for(int i=0;i< p->femp->esp->adp->NumOfEleY ;i++){
			for(int j=0;j< p->femp->esp->adp->NumOfEleY ;j++){
				double tmp_sigma = p->SigmaMap_Assume[p->femp->esp->adp->ConvertElement(j,i)];
				std::cout << tmp_sigma << ",";
			}
			std::cout << std::endl;
		}
		//Do Inverse
		std::cout << "[DoInverse()]" << std::endl;
		p->DoInverse();
	}

	return 0;	

}

int main(int argc,char *argv[]){

	InverseProblem *p = new InverseProblem(new Fem(new ESet(
		//(new AnalysisDomain())->Set(2,(int[]){8,8,0},(double[]){0.15,0.15,0.08},(int[]){0,0,1})
		(new AnalysisDomain())->Test2D_Exp1()
	)));

	p->Init();
	p->LoadPhantomZ();

	WINDOW *base_win;
	initscr();
	curs_set(0);
	start_color();
	init_pair(0, COLOR_WHITE, COLOR_BLACK);
	init_pair(1, COLOR_WHITE, COLOR_BLUE);
	init_pair(2, COLOR_WHITE, COLOR_CYAN);
	init_pair(3, COLOR_WHITE, COLOR_GREEN);
	init_pair(4, COLOR_WHITE,COLOR_MAGENTA);
	init_pair(5, COLOR_WHITE,COLOR_RED);
	init_pair(6, COLOR_WHITE,COLOR_WHITE);
	init_pair(7, COLOR_WHITE,COLOR_YELLOW);
	base_win = newwin(25,80,0,0);
	double color_max;
	double color_min;
	double color_step;
	for(int cnt=0;;cnt++){

		//Define Color Map
		color_max = p->SigmaMap_Assume[0];
		color_min = p->SigmaMap_Assume[0];
		for(int i=0;i<(p->femp->SumOfElement);i++){
			if(p->SigmaMap_Assume[i] > color_max) color_max = p->SigmaMap_Assume[i];
			if(p->SigmaMap_Assume[i] < color_min) color_min = p->SigmaMap_Assume[i];
		}
		color_step = (color_max-color_min) / 6.0;

		//TUI(main)
		werase(base_win);
		wattrset(base_win, COLOR_PAIR(0) | WA_BOLD);
		wclrscr(base_win);
		box(base_win, 0, 0);
		wCenterTitle(base_win, "SigmaMap");
		for(int i=0;i< p->femp->esp->adp->NumOfEleY ;i++){
			for(int j=0;j< p->femp->esp->adp->NumOfEleY ;j++){
				double tmp_sigma = p->SigmaMap_Assume[p->femp->esp->adp->ConvertElement(j,i)];
				//Calc Color
				int tmp_color;
				if(color_step==0.0) tmp_color = 0;
				else tmp_color = (int)((tmp_sigma-color_min)/color_step);
				//Show
				std::stringstream tmp;
				tmp << std::setw(3) << std::setprecision(2) << tmp_sigma;
				wattrset(base_win, COLOR_PAIR(tmp_color));
				mvwaddstr(base_win, i+1 ,j*3 + 1, tmp.str().c_str() );
			}
		}
		wrefresh(base_win);

		p->DoInverse();

	}
	endwin();
	
	return 0;
}
