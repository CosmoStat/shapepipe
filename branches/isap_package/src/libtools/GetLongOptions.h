/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: A. Woiselle
**
**    Date:  June, 11th 2008
**    
**    File:  GetLongOptions.h
**
**
*******************************************************************************
**
**    DESCRIPTION  Long Option input management
**    ----------- 
**                 
**    Options start with '--'             
**		ex : --verbose --value 10.5
**                 
**    Usage in the main : 
**		+ add at the begining of main.cc : 
**	      	#include "GetLongOptions.h"
**	      	char** short_argv;
**	      	int short_argc;
**	      	
**		+ add theese lines to replace "filtinit(argc, argv);" :
**	      	map<string, string> opts;
**	      	GetLongOptions(argv, argc, opts);
**	      	filtinit(short_argc, short_argv);
**
**	    + examples for each argument :
**	      		Bool 			GetLongOption(opts,"verbose",Verbose);
**				float/bool/...  GetLongOption(opts,"SigmaNoise",SigmaNoise);
**				other ex 		bool tmp; GetLongOption(opts,"hard",tmp) && (P.FilterType = FT_HARD);
**
******************************************************************************/
#ifndef _GETLONGOPTIONS_H_
#define _GETLONGOPTIONS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <map>
#include "GlobalInc.h"

using namespace std;

char** short_argv;
int short_argc;
inline bool is_digit(char a)
{
	int b=(int)a;
	return ((b<58) && (b>47));
}

int GetLongOptions(char** _argv, int _argc, map<string,string>& options)
{
	options.clear();

	string name; 
	
	int narg=_argc;
	
	string cont;
	string prefix("--");
	string shortprefix("-");
	size_t found,needval;
	
	int cnt=_argc;
	
	// Read the long arguments
	for(int i=1;i<narg;i++)
	{
		cont = string(_argv[i]);
		found = cont.find(prefix);
		if (found!=string::npos)
		{
			string next=string(_argv[i+1]);
			needval = next.find(shortprefix);
			// ( number[0] or -number[1] ) and not the last one
			if( i<narg-1 && ( is_digit(_argv[i+1][0]) || (needval==0 && is_digit(_argv[i+1][1])) ) ) 
			{
				//cerr<<"needval:"<<cont<<"->"<<next<<endl;
				i++;
				options[cont] = next;
				cnt-=2;
			}
			else // bool
			{
				options[cont] = string("1");
				cnt-=1;
			}
		}
	}
	short_argc=cnt;
	short_argv = new char*[short_argc];
	for(int i=0;i<short_argc;i++)
		short_argv[i] = new char[64];

	cnt=0;
	// keep the short ones
	for(int i=0;i<narg;i++)
	{
		cont = string(_argv[i]);
		found = cont.find(prefix);
		if (found!=string::npos)
		{
			string next=string(_argv[i+1]);
			needval = next.find(shortprefix);
			if(  i<narg-1 && (is_digit(_argv[i+1][0]) || (needval==0 && is_digit(_argv[i+1][1]))) )
				i++;
		}
		else
		{
			strcpy(short_argv[cnt],_argv[i]);
			cnt++;
		}
	}
	for(int i=cnt;i<short_argc;i++)
		strcpy(short_argv[i],_argv[_argc-short_argc+i]); // filein, mask, fileout ...

	return 0;
}


	
	
	
template <class T> bool GetLongOption(map<string,string>& opts, const char* optname, T &toto)
{
	char opt[128];
	map<string,string>::iterator it;
	strcpy(opt,"--");
	strcat(opt,optname);
	it = opts.find(opt);
	if(it!=opts.end()){ istringstream ss(it->second); ss>>toto; return true;}
	return false;
};

bool GetLongOption(map<string,string>& opts, const char* optname, bool &toto)
{
	char opt[128];
	map<string,string>::iterator it;
	strcpy(opt,"--");
	strcat(opt,optname);
	it = opts.find(opt);
	if(it!=opts.end()){toto=true; return true;}
	return false;
};

bool GetLongOption(map<string,string>& opts, const char* optname, Bool &toto)
{
	char opt[128];
	map<string,string>::iterator it;
	strcpy(opt,"--");
	strcat(opt,optname);
	it = opts.find(opt);
	if(it!=opts.end()){ bool plop; istringstream ss(it->second); ss>>plop; toto = (Bool)plop; return true;}
	return false;
};


	
#endif
