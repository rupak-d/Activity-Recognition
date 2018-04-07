#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <dirent.h>

using namespace std;

void InputFile();

vector<string>file;
#define PI 3.141592

int M = 5;    //Number of bins for Theta.
int N = 5;    //Number of bins for Distance.


int main()
{
  InputFile();
	int count_files;
	int tot = 1;                                        // **************1 for TRAIN & 2 for TEST******************
	if (tot == 1)
	{
		count_files = 192;
	}
	else 
	{
		count_files = 128;
	}

  for (int count = 0; count < count_files; count++)            
  {
	ifstream inFile;
	inFile.open(file[count].c_str());
	cout<<"Reading File: "<<file[count]<<endl;             // Shows file being processed.

	if (!inFile.is_open())
	{
	  cout<<"Unable to open file"<<endl;                   // Shows Error message if file cannot be opened.
	}

	int line = 0;
	string s;
	while(getline(inFile, s))
	{	
	 line++;
	}	
	
	float values[line][5];                                // Size of array required.

	inFile.clear();
	inFile.seekg(0);
	string str;
	int n = 0;
	int m = 0;
  	float num = 0;

	while(getline(inFile, str))
	{	
	  if (m>4)
	  { 
	  n = n+1;
	  m = 0;   
	  }	
	  	   stringstream stream(str);
			while (stream >> num)
			{
			  values[n][m] = num;
			 m = m+1;
			}	
		//n++;
	}

int frames = (n+1)/20;
inFile.close();	                                     // Close Input file.

	
//----------------------------------------------------------------------------------------------


//Calculating Joint Distances...
float dist[frames][5];
for (int fd=0; fd<frames; fd++)
{ 
  int r0 = (2 + 20*(fd))-1;
  int r1 = (4 + 20*(fd))-1;
  int r2 = (8 + 20*(fd))-1;
  int r3 = (12 + 20*(fd))-1;
  int r4 = (16 + 20*(fd))-1;
  int r5 = (20 + 20*(fd))-1;

	dist[fd][0] = pow((pow((values[r1][2]-values[r0][2]),2) + pow((values[r1][3]-values[r0][3]),2) + pow((values[r1][4]-values[r0][4]),2)),0.5);

	dist[fd][1] = pow((pow((values[r2][2]-values[r0][2]),2) + pow((values[r2][3]-values[r0][3]),2) + pow((values[r2][4]-values[r0][4]),2)),0.5);

	dist[fd][2] = pow((pow((values[r3][2]-values[r0][2]),2) + pow((values[r3][3]-values[r0][3]),2) + pow((values[r3][4]-values[r0][4]),2)),0.5);

	dist[fd][3] = pow((pow((values[r4][2]-values[r0][2]),2) + pow((values[r4][3]-values[r0][3]),2) + pow((values[r4][4]-values[r0][4]),2)),0.5);

	dist[fd][4] = pow((pow((values[r5][2]-values[r0][2]),2) + pow((values[r5][3]-values[r0][3]),2) + pow((values[r5][4]-values[r0][4]),2)),0.5);

   //cout<<dist[fd][0]<<" "<<dist[fd][1]<<" "<<dist[fd][2]<<" "<<dist[fd][3]<<" "<<dist[fd][4]<<endl;
   //cout<<fd<<endl;
}



//----------------------------------------------------------------------------------------------


//Calculating Angles...
float theta[frames][5];
float a = 0;
float b = 0;
float c = 0;
float cos_t = 0;

for (int ft=0; ft<frames; ft++)
{ 
 int r0 = (2 + 20*(ft))-1;
 int r1 = (4 + 20*(ft))-1;
 int r2 = (8 + 20*(ft))-1;
 int r3 = (12 + 20*(ft))-1;
 int r4 = (16 + 20*(ft))-1;
 int r5 = (20 + 20*(ft))-1;

  a = pow((values[r1][2]-values[r0][2]),2) + pow((values[r1][3]-values[r0][3]),2);
  a = pow(a,0.5);
  b = pow((values[r3][2]-values[r0][2]),2) + pow((values[r3][3]-values[r0][3]),2);
  b = pow(b,0.5);
  c = pow((values[r1][2]-values[r3][2]),2) + pow((values[r1][3]-values[r3][3]),2);
  c = pow(c,0.5);
  cos_t = (pow(a,2)+pow(b,2)-pow(c,2))/(2*a*b);
  if (cos_t < -1.0) cos_t = -1.0 ;
  if (cos_t > 1.0) cos_t = 1.0 ;
  theta[ft][0] = acos(cos_t)*180/PI;

  a = pow((values[r1][2]-values[r0][2]),2) + pow((values[r1][3]-values[r0][3]),2);
  a = pow(a,0.5);
  b = pow((values[r2][2]-values[r0][2]),2) + pow((values[r2][3]-values[r0][3]),2);
  b = pow(b,0.5);
  c = pow((values[r1][2]-values[r2][2]),2) + pow((values[r1][3]-values[r2][3]),2);
  c = pow(c,0.5);
  cos_t = (pow(a,2)+pow(b,2)-pow(c,2))/(2*a*b);
  if (cos_t < -1.0) cos_t = -1.0 ;
  if (cos_t > 1.0) cos_t = 1.0 ;
  if ((isnan(acos(cos_t)*180/PI))>0)
  theta[ft][1] = theta[ft-1][1];
  else theta[ft][1] = acos(cos_t)*180/PI;

  a = pow((values[r2][2]-values[r0][2]),2) + pow((values[r2][3]-values[r0][3]),2);
  a = pow(a,0.5);
  b = pow((values[r4][2]-values[r0][2]),2) + pow((values[r4][3]-values[r0][3]),2);
  b = pow(b,0.5);
  c = pow((values[r2][2]-values[r4][2]),2) + pow((values[r2][3]-values[r4][3]),2);
  c = pow(c,0.5);
  cos_t = (pow(a,2)+pow(b,2)-pow(c,2))/(2*a*b);
  if (cos_t < -1.0) cos_t = -1.0 ;
  if (cos_t > 1.0) cos_t = 1.0 ;
  if ((isnan(acos(cos_t)*180/PI))>0)
  theta[ft][2] = theta[ft-1][2];
  else theta[ft][2] = acos(cos_t)*180/PI;

  a = pow((values[r5][2]-values[r0][2]),2) + pow((values[r5][3]-values[r0][3]),2);
  a = pow(a,0.5);
  b = pow((values[r4][2]-values[r0][2]),2) + pow((values[r4][3]-values[r0][3]),2);
  b = pow(b,0.5);
  c = pow((values[r5][2]-values[r4][2]),2) + pow((values[r5][3]-values[r4][3]),2);
  c = pow(c,0.5);
  cos_t = (pow(a,2)+pow(b,2)-pow(c,2))/(2*a*b);
  if (cos_t < -1.0) cos_t = -1.0 ;
  if (cos_t > 1.0) cos_t = 1.0 ;
  theta[ft][3] = acos(cos_t)*180/PI;

  a = pow((values[r5][2]-values[r0][2]),2) + pow((values[r5][3]-values[r0][3]),2);
  a = pow(a,0.5);
  b = pow((values[r3][2]-values[r0][2]),2) + pow((values[r3][3]-values[r0][3]),2);
  b = pow(b,0.5);
  c = pow((values[r5][2]-values[r3][2]),2) + pow((values[r5][3]-values[r3][3]),2);
  c = pow(c,0.5);
  cos_t = (pow(a,2)+pow(b,2)-pow(c,2))/(2*a*b);
  if (cos_t < -1.0) cos_t = -1.0 ;
  if (cos_t > 1.0) cos_t = 1.0 ;
  theta[ft][4] = acos(cos_t)*180/PI;
//cout<<theta[ft][0]<<" "<<theta[ft][1]<<" "<<theta[ft][2]<<" "<<theta[ft][3]<<" "<<theta[ft][4]<<endl;
//cout<<ft<<endl;
}


//----------------------------------------------------------------------------------------------	


//Histograms...
float d_max[5] = {0};
float t_max[5] = {0};
float d_min[5] = {0};
float t_min[5] = {0};
float d_range[5] = {0};
float t_range[5] = {0};
float d_binsize[5] = {0};
float t_binsize[5] = {0};

for (int w=0; w < 5; w++)
{
    for (int q = 0; q < frames; q++)
    {
      if (dist[q][w] > d_max[w])
      {
      d_max[w] = dist[q][w];
      }
      if (theta[q][w] > t_max[w])
      {
      t_max[w] = theta[q][w];
      }
	  if (dist[q][w] < d_min[w])
      {
      d_min[w] = dist[q][w];
      }
	  if (theta[q][w] < t_min[w])
      {
      t_min[w] = theta[q][w];
      }
    }
   d_range[w] = d_max[w] - d_min[w];
   t_range[w] = t_max[w] - t_min[w];
   d_binsize[w] = d_range[w]/N;
   t_binsize[w] = t_range[w]/M;

}

float ref_dist[frames][5];
float ref_theta[frames][5];


	for (int w = 0; w < 5; w++)
    {
	  for (int q = 0; q < frames; q++)
      {
	  ref_dist[q][w] = dist[q][w] - d_min[w];
	  ref_theta[q][w] = theta[q][w] - t_min[w];
	  }
	}


vector<float> histogram((5*(M+N)),0);

	for (int q = 0; q < frames; q++)
    {  
	 int bin_number[10] = {0};
	 for (int w = 0; w < 5; w++)
	 {
	  bin_number[w] = int (ceil(ref_dist[q][w]/d_binsize[w]));
	  bin_number[w+N] = int (ceil(ref_theta[q][w]/t_binsize[w]));
		if(bin_number[w] == 0)
		{
          bin_number[w] = 1;
        }
        if(bin_number[w+N] == 0)
		{
          bin_number[w+N] = 1;
        }
	 }
	  for(int w=0;w<5;w++)
	  {
	   histogram[(w*N)+((bin_number[w])-1)] = histogram[(w*N)+((bin_number[w])-1)] + 1;
	   histogram[((5*N)+(w*5))+((bin_number[5+w])-1)] = histogram[((5*N)+(w*5))+((bin_number[5+w])-1)] + 1;
	  }
    }


//------------------------------------------------------------------------------------------------


//Output File...
	ofstream outFile;
	int se = 0;
	outFile.open("rad_d1",ofstream::app);                   // Uncomment if Train files are being processed.
	//outFile.open("rad_d1.t",ofstream::app);               // Uncomment if Test files are being processed.

	if(tot==1)
	{
	se = 12;
	}
	else
	{
	se = 8;
	}
	for (int i=0;i<5*(M+N);i++)
	{
		if (i == 0)
		{
		 outFile << (count/se)+1 << " ";                        
		}               
     	outFile << i+1 << ":" << histogram[i]/frames<<" ";
	}
	outFile << "\n";
outFile.close();
}

return 0;	
}


//-------------------------------------------------------------------------------------------------


//File Names..
void InputFile()
{
ofstream files("filenames.txt");
	
//************TRAINING SET FILES****************************************************************

  for (int a=1; a<17; a++)
  {
	for (int s=1; s<7; s++)
	{
	  for (int e=1; e<3; e++)
	  {
		if(a<10)
		files << "../dataset_full/Train/a0"<<a<<"_s0"<<s<<"_e0"<<e<<"_skeleton_proj.txt"<<endl;
		else 
		files << "../dataset_full/Train/a"<<a<<"_s0"<<s<<"_e0"<<e<<"_skeleton_proj.txt"<<endl;
	  }
	}
  }
//*********************************************************************************************

//************TEST SET FILES*******************************************************************

  /*for (int a=1; a<17; a++)
  {
	for (int s=7; s<11; s++)
	{
	  for (int e=1; e<3; e++)
	  {
		if(a>9 && s<10)
		files << "../dataset_full/Test/a"<<a<<"_s0"<<s<<"_e0"<<e<<"_skeleton_proj.txt"<<endl;
		else if(a<10 && s>9)
		files << "../dataset_full/Test/a0"<<a<<"_s"<<s<<"_e0"<<e<<"_skeleton_proj.txt"<<endl;
		else if(a>9 && s>9)
		files << "../dataset_full/Test/a"<<a<<"_s"<<s<<"_e0"<<e<<"_skeleton_proj.txt"<<endl;
		else 
		files << "../dataset_full/Test/a0"<<a<<"_s0"<<s<<"_e0"<<e<<"_skeleton_proj.txt"<<endl;
	  }
	}
  }*/
//*********************************************************************************************

 files.close();

	ifstream readFile;
	readFile.open("filenames.txt");                 // Creates file with names of all Input files.
	string str;
	string temp;
	while(getline(readFile, str))
	{	
	  	stringstream stream(str);
			while (stream >> temp)
			{
			 file.push_back(temp);
			}	
	}
 readFile.close();                          // Close File.
}
