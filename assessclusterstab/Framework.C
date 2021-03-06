#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readClusters(const char* aFName)
{
	clusterings=0;
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		readCluster(buffer);
	}
	inFile.close();
	return 0;
}

int 
Framework::showAdjacencyMap(const char* aFName)
{
	ofstream oFile(aFName);
	oFile <<"ORF";
	for(map<string,map<string,double>*>::iterator vIter=adjacencymatrix.begin();vIter!=adjacencymatrix.end();vIter++)
	{
		oFile << "\t" << vIter->first;
	}
	oFile <<endl;
	double maxclusterings=0;
	for(map<string,map<string,double>*>::iterator vIter=adjacencymatrix.begin();vIter!=adjacencymatrix.end();vIter++)
	{
		map<string,double>* edgeStrength=vIter->second;
		for(map<string,map<string,double>*>::iterator uIter=adjacencymatrix.begin();uIter!=adjacencymatrix.end();uIter++)
		{
			double val=0;
			if(edgeStrength->find(uIter->first)!=edgeStrength->end())
			{
				val=(*edgeStrength)[uIter->first];
				if(val>maxclusterings)
				{
					maxclusterings=val;
				}
			}
		}
	}
	for(map<string,map<string,double>*>::iterator vIter=adjacencymatrix.begin();vIter!=adjacencymatrix.end();vIter++)
	{
		oFile <<vIter->first;
		map<string,double>* edgeStrength=vIter->second;
		double total=(*edgeStrength)[vIter->first];
		for(map<string,map<string,double>*>::iterator uIter=adjacencymatrix.begin();uIter!=adjacencymatrix.end();uIter++)
		{
			double val=0;
			if(edgeStrength->find(uIter->first)!=edgeStrength->end())
			{
				//val=(*edgeStrength)[uIter->first]/maxclusterings;
				val=(*edgeStrength)[uIter->first]/total;
			}
			oFile <<"\t" << val;
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}


int
Framework::readCluster(const char* aFName)
{
	clusterings++;
	map<int,map<string,int>*> clusterSet; 
	char buffer[1024];
	ifstream inFile(aFName);
	
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string geneName;
		int clusterid=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else if(tokCnt==1)
			{
				clusterid=atoi(tok);
			}		
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* cluster=NULL;
		if(clusterSet.find(clusterid)==clusterSet.end())
		{
			cluster=new map<string,int>;
			clusterSet[clusterid]=cluster;
		}
		else
		{
			cluster=clusterSet[clusterid];
		}
		(*cluster)[geneName]=0;
	}
	for(map<int,map<string,int>*>::iterator cIter=clusterSet.begin();cIter!=clusterSet.end();cIter++)
	{
		map<string,int>* cluster=cIter->second;
		for(map<string,int>::iterator gIter=cluster->begin();gIter!=cluster->end();gIter++)
		{
			map<string,double>* generow_g=NULL;
			if(adjacencymatrix.find(gIter->first)==adjacencymatrix.end())
			{
				generow_g=new map<string,double>;
				adjacencymatrix[gIter->first]=generow_g;
				(*generow_g)[gIter->first]=1;
			}
			else
			{
				generow_g=adjacencymatrix[gIter->first];
				(*generow_g)[gIter->first]=(*generow_g)[gIter->first]+1;
			}
			map<string,int>::iterator hIter=gIter;
			hIter++;
			for(;hIter!=cluster->end();hIter++)
			{
				map<string,double>* generow_h=NULL;
				if(adjacencymatrix.find(hIter->first)==adjacencymatrix.end())
				{
					generow_h=new map<string,double>;
					adjacencymatrix[hIter->first]=generow_h;
					(*generow_h)[hIter->first]=1;
				}
				else
				{
					generow_h=adjacencymatrix[hIter->first];
				}
				if(generow_g->find(hIter->first)==generow_g->end())
				{
					(*generow_g)[hIter->first]=1;
				}
				else
				{
					(*generow_g)[hIter->first]=(*generow_g)[hIter->first]+1;
				}
				if(generow_h->find(gIter->first)==generow_h->end())
				{
					(*generow_h)[gIter->first]=1;
				}
				else
				{
					(*generow_h)[gIter->first]=(*generow_h)[gIter->first]+1;
				}
			}
		}
	}
	inFile.close();
	for(map<int,map<string,int>*>::iterator cIter=clusterSet.begin();cIter!=clusterSet.end();cIter++)
	{
		cIter->second->clear();
		delete cIter->second;
	}
	clusterSet.clear();
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=3)
	{
		cout <<"assessClusterStab clusterings outputfile" << endl;
		return 0;
	}
	Framework fw;
	fw.readClusters(argv[1]);
	fw.showAdjacencyMap(argv[2]);
	return 0;
}
