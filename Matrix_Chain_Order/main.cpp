//
//  main.cpp
//  Matrix_Chain_Order
//
//  Created by LI YANZHE on 15/07/2017.
//  Copyright © 2017 Yanzhe Lee. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <random>

void Show_Index_Page()
{
	puts("\n--------------------------------------------------------------------------------");
	puts("|                                                                              |");
	puts("|         Copyright (C) 2016 Yanzhe Lee. All rights reserved.                  |");
	puts("|                                                                              |");
	puts("|                       Harbin Institute of Technology                         |");
	puts("|                                                                              |");
	puts("|         License GPLv3+: GNU GPL version 3 or later                           |");
	puts("|                                                                              |");
	puts("|         This is free software: you are free to change and redistribute it    |");
	puts("|                                                                              |");
	puts("|         Email: lee.yanzhe@yanzhe.org           Version 1.0.0                 |");
	puts("|                                                                              |");
	puts("--------------------------------------------------------------------------------");
	puts("          Please maximize your window to get a better display effect            ");
	puts("--------------------------------------------------------------------------------");
	puts("|                                                                              |");
	puts("|       Matrix Multiplication Chain Optimal Parenthesis Solution Finder        |");
	puts("|                                                                              |");
	puts("--------------------------------------------------------------------------------");
}

class Matrix_Chain
{
public:
	Matrix_Chain(std::vector <int> p)
	{
		this->originChainDesc=p;
		this->chainLength=p.size()-1;
		this->minCostMap=Create_Matrix<int>(chainLength, chainLength);
		this->divisionPointMap=Create_Matrix<int>(chainLength, chainLength);
		this->operationFinish=false;
	}
	template <class T>
	static std::vector<std::vector<T>> Create_Matrix(unsigned long m,unsigned long n)
	{
		std::vector<std::vector<T>> matrix(m);
		for (size_t i=0; i<m; ++i) {
			matrix[i].resize(n,0);
		}
		return matrix;
	}
	
	template <class T>
	static void Show_Matrix(std::vector<std::vector<T>> matrix)
	{
		using namespace std;
		size_t m=matrix.size();
		size_t n=(m==0)?0:matrix[0].size();
		for (size_t i=0; i<m; ++i) {
			for (size_t j=0; j<n; ++j) {
				printf("%5d\t",matrix.at(i).at(j));
				//				cout<<matrix.at(i).at(j)<<"\t";
			}
			cout<<endl;
		}
	}
	unsigned long cost()
	{
		return cost_between(0,chainLength-1);
	}
	std::vector<std::vector<int>> cost_map()
	{
		check();
		return minCostMap;
	}
	
	int cost_between(int i,unsigned long j)
	{
		check();
		return minCostMap.at(i).at(j);
	}
	std::string parens()
	{
		return parens(0,chainLength-1);
	}
	std::string parens(unsigned long i,unsigned long j)
	{
		check();
		using namespace std;
		static string solution="";
		if (i==j) {
			solution.append("A");
		}else{
			solution.append("(");
			parens(i, divisionPointMap.at(i).at(j));
			parens(divisionPointMap.at(i).at(j)+1, j);
			solution.append(")");
		}
		return solution;
	}
	void print_chain()
	{
		using namespace std;
		size_t i;
		for (i=0; i<chainLength-2; ++i) {
			printf("%dx%d,",originChainDesc.at(i),originChainDesc.at(i+1));
		}
		printf("%dx%d",originChainDesc.at(i),originChainDesc.at(i+1));
		cout<<endl;
	}
	
protected:
	std::vector <int> originChainDesc;
	std::vector<std::vector<int>> minCostMap;
	std::vector<std::vector<int>> divisionPointMap;
	unsigned long chainLength;
	bool operationFinish;
	void Matrix_Chain_Order()
	{
		using namespace std;
		//		for (int i=0; i<chainLength; ++i) {
		//			minCostMap[i][i]=0;
		//		}
		for (int len=2; len<=chainLength; ++len) {
			for (int i=0; i<chainLength-len+1; ++i) {
				int j=i+len-1;
				minCostMap[i][j]=INT_MAX;
				for (int k=i; k<j; ++k) {
					int c=minCostMap[i][k]+minCostMap[k+1][j]+originChainDesc[i]*originChainDesc[k+1]*originChainDesc[j+1];
					if(c<minCostMap[i][j]){
						minCostMap[i][j]=c;
						divisionPointMap[i][j]=k;
					}
				}
			}
		}
		operationFinish=true;
	}
	inline void check()
	{
		if (!operationFinish)
			Matrix_Chain_Order();
	}
};


int main(int argc, const char * argv[]) {
	Show_Index_Page();
	using namespace std;
	char mode;
	bool testFlag=false;
	vector <int> scaleDescr;
	cout<<"Choose mode(I:input,T:mass_test):";
	if(
	   scanf("%c",&mode)!=1
	   ||
	   (mode!='I'&&mode!='i'&&mode!='T'&&mode!='t')
	   ){
		cout<<"Mode Error"<<endl;
		return 1;
	}
	else{
		switch (mode) {
			case 'T':
			case 't':testFlag=true;break;
			default:break;
		}
	}
	cout<<"\nPlease input Matrix-Chain length:";
	int n=0;
	cin>>n;
	if (!n) {
		cout<<"Length Error"<<endl;
		return 1;
	}
	if(!testFlag)cout<<"\nPlease input Scale Description:\n";
	std::default_random_engine gener;
	std::uniform_int_distribution<int> dis(1,1000);
	for (size_t i=0; i<n+1; ++i) {
		int tp=0;
		if(testFlag)tp=dis(gener);
		else cin>>tp;
		scaleDescr.push_back(tp);
	}
	Matrix_Chain Chain_Instance(scaleDescr);
	double startop = clock();
	cout<<"\n=============================    Origin Chain     ==============================\n\n";
	Chain_Instance.print_chain();
	//	Matrix_Chain::Show_Matrix(Chain_Instance.cost_map());
	cout<<"\n=============================     Total Cost      ==============================\n\n";
	printf("%lu Scalar Multiplication Operations\n",Chain_Instance.cost());
	cout<<"\n============================= Optimal Parenthesis ==============================\n\n";
	cout<<Chain_Instance.parens()<<endl;
	cout<<"\n=============================   Operation Time    ==============================\n\n";
	double endop = clock();
	printf("%lf s\n\n",(endop-startop)/CLOCKS_PER_SEC);
	
	return 0;
}
