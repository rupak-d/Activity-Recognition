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
float calc_dist(float x0, float y0, float x1, float y1);     // Calculates distance.
float comp_angle(float x0, float y0, float x1, float y1);    // Computes angle.

vector<string>file;
#define PI 3.141592
int number_of_quantised_theta = 10;                      // Divided into 8 bins.
float theta_range = 360/number_of_quantised_theta;      // Bin range of theta.



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
  for (int count = 0; count < test_files; count++)       
  {
	ifstream inFile;
	inFile.open(file[count].c_str());
	cout<<"Reading File: "<<file[count]<<endl;        // Shows file being processed.

	if (!inFile.is_open())
	{
	  cout<<"Unable to open file"<<endl;              // Shows Error message if file cannot be opened.
	}

	int line = 0;
	string s;
	while(getline(inFile, s))
	{	
	 line++;
	}	
	
	float values[line][5];                           // Size of array required.

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
inFile.close();	                                         // Close Input file.

	
//------------------------------------------------------------------------------------------------------


//Initialize Vectors...
vector<vector<vector<float> > > p_xy (frames,vector<vector<float> >(20,vector <float>(2,0)));
vector<vector<vector<float> > > p_yz (frames,vector<vector<float> >(20,vector <float>(2,0)));
vector<vector<vector<float> > > p_xz (frames,vector<vector<float> >(20,vector <float>(2,0)));

int number_of_instances = frames-1;
int level_2_1_instances = number_of_instances/2;
int level_2_2_instances = number_of_instances-level_2_1_instances; 
int level_2_1_1_instances = level_2_1_instances/2;
int level_2_1_2_instances = level_2_1_instances - level_2_1_1_instances;
int level_2_2_1_instances = level_2_2_instances/2;
int level_2_2_2_instances = level_2_2_instances - level_2_2_1_instances;

vector<vector<float> > d_xy(20, vector<float>(number_of_instances,0));
vector<vector<float> > d_yz(20, vector<float>(number_of_instances,0));
vector<vector<float> > d_xz(20, vector<float>(number_of_instances,0));

vector<vector<float> > t_xy(20, vector<float>(number_of_instances,0));
vector<vector<float> > t_yz(20, vector<float>(number_of_instances,0));
vector<vector<float> > t_xz(20, vector<float>(number_of_instances,0));

vector<vector<float> > q_xy(20, vector<float>(number_of_instances,0));
vector<vector<float> > q_yz(20, vector<float>(number_of_instances,0));
vector<vector<float> > q_xz(20, vector<float>(number_of_instances,0));

vector<vector<float> > hist_xy(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz(20,vector <float>(number_of_quantised_theta,0));

vector<vector<float> > hist_xy_level_2_1(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz_level_2_1(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz_level_2_1(20,vector <float>(number_of_quantised_theta,0));
		
vector<vector<float> > hist_xy_level_2_2(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz_level_2_2(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz_level_2_2(20,vector <float>(number_of_quantised_theta,0));
		
vector<vector<float> > hist_xy_level_2_1_1(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz_level_2_1_1(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz_level_2_1_1(20,vector <float>(number_of_quantised_theta,0));
		
vector<vector<float> > hist_xy_level_2_1_2(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz_level_2_1_2(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz_level_2_1_2(20,vector <float>(number_of_quantised_theta,0));

vector<vector<float> > hist_xy_level_2_2_1(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz_level_2_2_1(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz_level_2_2_1(20,vector <float>(number_of_quantised_theta,0));

vector<vector<float> > hist_xy_level_2_2_2(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_yz_level_2_2_2(20,vector <float>(number_of_quantised_theta,0));
vector<vector<float> > hist_xz_level_2_2_2(20,vector <float>(number_of_quantised_theta,0));

vector<vector<float> > hist_3D(20,vector <float>(21*number_of_quantised_theta,0));
vector<float> hist_3D_final((20*21*number_of_quantised_theta),0);

//------------------------------------------------------------------------------------------------------


//Filling in values...
for(int u=0; u<frames; u++)
{ 
	for(int i=0;i<20;i++)
	{
		for(int n=0;n<2;n++)
		{
			p_xy[u][i][n]=values[(20*u)+i][2+n];    
			p_yz[u][i][n]=values[(20*u)+i][3+n];    
			p_xz[u][i][n]=values[(20*u)+i][2+(2*n)];    
		}	
	}
}
    
	for(int m=0;m<number_of_instances; m++)
	{
		for(int i=0;i<20;i++)
		{
			d_xy[i][m]=calc_dist(p_xy[m][i][0],p_xy[m][i][1],p_xy[m+1][i][0],p_xy[m+1][i][1]);
			d_yz[i][m]=calc_dist(p_yz[m][i][0],p_yz[m][i][1],p_yz[m+1][i][0],p_yz[m+1][i][1]);
			d_xz[i][m]=calc_dist(p_xz[m][i][0],p_xz[m][i][1],p_xz[m+1][i][0],p_xz[m+1][i][1]);
			
			t_xy[i][m]=floor(comp_angle(p_xy[m][i][0],p_xy[m][i][1],p_xy[m+1][i][0],p_xy[m+1][i][1])*100)/100;
			t_yz[i][m]=floor(comp_angle(p_yz[m][i][0],p_yz[m][i][1],p_yz[m+1][i][0],p_yz[m+1][i][1])*100)/100;
			t_xz[i][m]=floor(comp_angle(p_xz[m][i][0],p_xz[m][i][1],p_xz[m+1][i][0],p_xz[m+1][i][1])*100)/100;		
		}	
	}


for(int m=0;m<number_of_instances;m++)
{
  for(int i=0;i<20;i++)
	{
			if(t_xy[i][m]<0){
				t_xy[i][m]=(t_xy[i][m])+360;
			}
			if(t_yz[i][m]<0){
				t_yz[i][m]=(t_yz[i][m])+360;
			}
			if(t_xz[i][m]<0){
				t_xz[i][m]=(t_xz[i][m])+360;
			}
	}
}

for(int m=0;m<number_of_instances;m++)
	{
		for(int i=0;i<20;i++)
		{
			q_xy[i][m]=ceil(t_xy[i][m]/theta_range);
			q_yz[i][m]=ceil(t_yz[i][m]/theta_range);
			q_xz[i][m]=ceil(t_xz[i][m]/theta_range);

			if(q_xy[i][m]==0){
				q_xy[i][m]=1;                           
			}
			if(q_yz[i][m]==0){
				q_yz[i][m]=1;
			}
			if(q_xz[i][m]==0){
				q_xz[i][m]=1;
			}
		}
	}


//---------------------------------------------------------------------------------------------------------------------------


// Making Histograms, Temporal Pyramids...
for(int m=0;m<number_of_instances;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy[i][q_xy[i][m]-1] = hist_xy[i][q_xy[i][m]-1] + d_xy[i][m];
			hist_yz[i][q_yz[i][m]-1] = hist_yz[i][q_xy[i][m]-1] + d_yz[i][m];
			hist_xz[i][q_xz[i][m]-1] = hist_xz[i][q_xy[i][m]-1] + d_xz[i][m];
		}
}


for(int m=0;m< level_2_1_instances ;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy_level_2_1[i][q_xy[i][m]-1]  =  hist_xy_level_2_1[i][q_xy[i][m]-1]  +  d_xy[i][m];
			hist_yz_level_2_1[i][q_yz[i][m]-1]  =  hist_yz_level_2_1[i][q_xy[i][m]-1]  +  d_yz[i][m];
			hist_xz_level_2_1[i][q_xz[i][m]-1]  =  hist_xz_level_2_1[i][q_xy[i][m]-1]  + d_xz[i][m];
		}
}
for(int m=0;m< level_2_2_instances ;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy_level_2_2[i][q_xy[i][m]-1] = hist_xy_level_2_2[i][q_xy[i][m]-1] + d_xy[i][m + level_2_1_instances];
			hist_yz_level_2_2[i][q_yz[i][m]-1] = hist_yz_level_2_2[i][q_yz[i][m]-1] + d_yz[i][m + level_2_1_instances];
			hist_xz_level_2_2[i][q_xz[i][m]-1] = hist_xz_level_2_2[i][q_xz[i][m]-1] + d_xz[i][m + level_2_1_instances];
		}
}
for(int m=0;m< level_2_1_1_instances ;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy_level_2_1_1[i][q_xy[i][m]-1] = hist_xy_level_2_1_1[i][q_xz[i][m]-1] + d_xy[i][m];
			hist_yz_level_2_1_1[i][q_yz[i][m]-1] = hist_yz_level_2_1_1[i][q_xz[i][m]-1] + d_yz[i][m];
			hist_xz_level_2_1_1[i][q_xz[i][m]-1] = hist_xz_level_2_1_1[i][q_xz[i][m]-1] + d_xz[i][m];
		}
}
for(int m=0;m< level_2_1_2_instances ;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy_level_2_1_2[i][q_xy[i][m]-1] = hist_xy_level_2_1_2[i][q_xy[i][m]-1] + d_xy[i][m + level_2_1_1_instances];
			hist_yz_level_2_1_2[i][q_yz[i][m]-1] = hist_yz_level_2_1_2[i][q_xy[i][m]-1] + d_yz[i][m + level_2_1_1_instances];
			hist_xz_level_2_1_2[i][q_xz[i][m]-1] = hist_xz_level_2_1_2[i][q_xy[i][m]-1] + d_xz[i][m + level_2_1_1_instances];
		}
}
for(int m=0;m< level_2_2_1_instances ;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy_level_2_2_1[i][q_xy[i][m]-1] = hist_xy_level_2_2_1[i][q_xy[i][m]-1] + d_xy[i][m + level_2_1_instances];
			hist_yz_level_2_2_1[i][q_yz[i][m]-1] = hist_yz_level_2_2_1[i][q_xy[i][m]-1] + d_yz[i][m + level_2_1_instances];
			hist_xz_level_2_2_1[i][q_xz[i][m]-1] = hist_xz_level_2_2_1[i][q_xy[i][m]-1] + d_xz[i][m + level_2_1_instances];
		}
}
for(int m=0;m< level_2_2_2_instances ;m++)
{
		for(int i=0;i<20;i++)
		{
			hist_xy_level_2_2_2[i][q_xy[i][m]-1] = hist_xy_level_2_2_2[i][q_xy[i][m]-1] + d_xy[i][m + level_2_1_instances + level_2_2_1_instances];
			hist_yz_level_2_2_2[i][q_yz[i][m]-1] = hist_yz_level_2_2_2[i][q_xy[i][m]-1] + d_yz[i][m + level_2_1_instances + level_2_2_1_instances];
			hist_xz_level_2_2_2[i][q_xz[i][m]-1] = hist_xz_level_2_2_2[i][q_xy[i][m]-1] + d_xz[i][m + level_2_1_instances + level_2_2_1_instances];
		}
}


for(int i=0;i<20;i++)
{
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][j] = hist_xy[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][j+number_of_quantised_theta] = hist_xy_level_2_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(2*number_of_quantised_theta) + j] = hist_xy_level_2_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(3*number_of_quantised_theta)+j] = hist_xy_level_2_1_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(4*number_of_quantised_theta)+j] = hist_xy_level_2_1_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(5*number_of_quantised_theta)+j] = hist_xy_level_2_2_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(6*number_of_quantised_theta)+j] = hist_xy_level_2_2_2[i][j];  
		}

		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(7*number_of_quantised_theta)+j] = hist_yz[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(8*number_of_quantised_theta)+j] = hist_yz_level_2_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(9*number_of_quantised_theta)+j] = hist_yz_level_2_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(10*number_of_quantised_theta)+j] = hist_yz_level_2_1_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(11*number_of_quantised_theta)+j] = hist_yz_level_2_1_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(12*number_of_quantised_theta)+j] = hist_yz_level_2_2_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(13*number_of_quantised_theta)+j] = hist_yz_level_2_2_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(14*number_of_quantised_theta)+j] = hist_xz[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(15*number_of_quantised_theta)+j] = hist_xz_level_2_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(16*number_of_quantised_theta)+j] = hist_xz_level_2_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(17*number_of_quantised_theta)+j] = hist_xz_level_2_1_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(18*number_of_quantised_theta)+j] = hist_xz_level_2_1_2[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(19*number_of_quantised_theta)+j] = hist_xz_level_2_2_1[i][j];  
		}
		for(int j=0;j < number_of_quantised_theta; j++){
				hist_3D[i][(20*number_of_quantised_theta)+j] = hist_xz_level_2_2_2[i][j];  
		}

}


for(int i=0;i<20;i++)
{
	for(int j=0;j<(21*number_of_quantised_theta);j++)
	{
		hist_3D_final[(i*(21*number_of_quantised_theta))+j] = hist_3D[i][j];
	}
} 


//----------------------------------------------------------------------------------------------------------------------


//Output File...
	ofstream outFile;
	//outFile.open("hod_d2",ofstream::app);                  // Uncomment if Train files are being processed.
	outFile.open("hod_d2.t",ofstream::app);              // Uncomment if Test files are being processed.
	int se = 0;
	if(tot==1)
	{
	se = 12;
	}
	else
	{
	se = 8;
	}
	for (int i=0;i<(21*20*number_of_quantised_theta);i++)
	{
	 if (i==0){
	 outFile << (count/se)+1 << " ";}                    
		if (i<((21*20*number_of_quantised_theta)-1))
		{
		outFile << i+1 << ":" << hist_3D_final[i];
		outFile << " ";
		}
		else
		{	
		outFile << i+1 << ":" << hist_3D_final[i];
		outFile << "\n";
		}
	}

outFile.close();

	}
	return 0;
	
}


//----------------------------------------------------------------------------------------------------------------


// Distance function.
float calc_dist(float x0, float y0, float x1, float y1)
{
		float x_sqr = (x0 - x1)*(x0 - x1);
		float y_sqr = (y0 - y1)*(y0 - y1);
		float dist_sqr = x_sqr+y_sqr;
		return sqrt(dist_sqr);
}
	
// Angle function.
float comp_angle(float x0, float y0, float x1, float y1){
			float angle= atan2((y1-y0),(x1-x0));
			angle =angle*180/PI;
			return (angle);
}


//------------------------------------------------------------------------------------------------------------


//File Names..
void InputFile()
{
ofstream files("filenames.txt");
	
//************TRAINING SET FILES*****************************************************************************

  /*for (int a=1; a<17; a++)
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
  }*/
//************************************************************************************************************

//*************TEST SET FILES*********************************************************************************

  for (int a=1; a<17; a++)
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
  }
//************************************************************************************************************

 files.close();

	ifstream readFile;
	readFile.open("filenames.txt");               // Creates file with names of all Input files.
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
 readFile.close();                                   // Close File.
}










